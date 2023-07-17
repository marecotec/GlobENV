##################################################################################################
# GlobENV - Data Validator - GlobENV_Validator.py
#
# Dr Andy Davies/Martyn Roberts - Bangor University, Wales - github: marecotec
#
# This python code requires ArcGIS 10.4+.
#
##################################################################################################
# Copyright 2018 Dr Andy Davies
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished
# to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
# NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##################################################################################################


# PURPOSE
# Script that reads in bottom of cast (global) CTD dataset (post 2013), produced by first running ODV_WOD_CTD_ExtractCoordsandDeepestSample.py,
# and uses it to validate global environmental seabed layers produced by Andrew Davies' tri-linear WOA interpolation approach (WOA is pre 2013 data).
# Works much faster when inputs are GRID files rather than ASCII (don't know why).
# MAKE SURE YOU ONLY PUT LAYERS YOU WANT TO VALIDATE IN 'TRI_PATH' - In this case, temp, sal, diso2.

# IMPORT PACKAGES
import arcpy
from arcpy.sa import *
import unicodedata
import csv
import numpy
import os
import matplotlib.pyplot as plt
from scipy.stats import linregress, stats
import pylab as lab
import time
import pandas as pd
from collections import namedtuple
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gc

arcpy.CheckOutExtension("Spatial")
pd.options.mode.chained_assignment = None


def SpatialErrorMap(depth_cut_off, wod_var_name):
    spatial_error_array = BuildDFArcFC(os.path.join(output_directory, "Shapefiles", "CTD_location_points_p_sj.shp"))
    spatial_error_array = spatial_error_array[spatial_error_array[wod_var_name] != -999999]
    spatial_error_array = spatial_error_array[spatial_error_array[wod_var_name] != -9999]

    spatial_error_array['dep_diff'] = spatial_error_array['Depth'] / spatial_error_array['bathy_dep']

    spatial_error_array_0 = spatial_error_array[spatial_error_array['dep_diff'] >= 1.0 - depth_cut_off]
    spatial_error_array_1 = spatial_error_array_0[spatial_error_array_0['dep_diff'] <= 1.0 + depth_cut_off]

    spatial_error_array_1['tri_diff'] = spatial_error_array_1[wod_var_name] - spatial_error_array_1['tri']

    spatial_error_array_2 = spatial_error_array_1.groupby('Id')['tri_diff'].agg(['count', 'mean']).reset_index()

    id_count_dict = dict(zip(spatial_error_array_2['Id'], spatial_error_array_2['count']))
    id_mean_dict = dict(zip(spatial_error_array_2['Id'], spatial_error_array_2['mean']))

    if not os.path.exists(os.path.join(output_directory, "Plots", "spatial_error.png")):

        inRas = arcpy.Raster(os.path.join(output_directory, "Rasters", "fishnet"))
        xmin = inRas.extent.XMin
        ymin = inRas.extent.YMin
        xmax = inRas.extent.XMax
        ymax = inRas.extent.YMax

        nrows = int((inRas.extent.YMax - inRas.extent.YMin) / inRas.meanCellHeight)
        ncols = int((inRas.extent.XMax - inRas.extent.XMin) / inRas.meanCellWidth)

        del inRas

        fig = plt.figure()
        fig.subplots_adjust(hspace=0.1, wspace=0.1)

        grd = arcpy.RasterToNumPyArray(os.path.join(output_directory, "Rasters", "fishnet"), arcpy.Point(xmin, ymin),
                                       ncols, nrows)
        grd2 = numpy.empty(grd.shape, dtype=float)
        grd2[:] = numpy.nan
        for k, v in id_count_dict.items():
            grd2[grd == k] = v

        ax = fig.add_subplot(1, 2, 1)
        ax.set_title("Count of casts (a)")
        im = ax.matshow(grd2, extent=[xmin, xmax, ymin, ymax])
        ax.set_xticks([])
        ax.set_yticks([])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        plt.colorbar(im, cax=cax)
        del grd, grd2

        grd = arcpy.RasterToNumPyArray(os.path.join(output_directory, "Rasters", "fishnet"), arcpy.Point(xmin, ymin),
                                       ncols, nrows)
        grd2 = numpy.empty(grd.shape, dtype=float)
        grd2[:] = numpy.nan
        for k, v in id_mean_dict.items():
            grd2[grd == k] = v

        ax = fig.add_subplot(1, 2, 2)
        ax.set_title("Mean difference (b)")
        im = ax.matshow(grd2, extent=[xmin, xmax, ymin, ymax])
        ax.set_xticks([])
        ax.set_yticks([])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        plt.colorbar(im, cax=cax)
        del grd, grd2

        fig.savefig(os.path.join(output_directory, "Plots", "spatial_error.png"), dpi=300, bbox_inches='tight')
        del fig, ax


def ArcpyMap_1(input_bathymetry, input_tri, output_directory):
    input_layer_prop = {}

    input_layer_prop["cell size x"] = arcpy.GetRasterProperties_management(input_tri,
                                                                           "CELLSIZEX").getOutput(0)
    input_layer_prop["cell size y"] = arcpy.GetRasterProperties_management(input_tri,
                                                                           "CELLSIZEY").getOutput(0)
    input_layer_prop["nrows"] = arcpy.GetRasterProperties_management(input_tri,
                                                                     "ROWCOUNT").getOutput(0)
    input_layer_prop["ncols"] = arcpy.GetRasterProperties_management(input_tri,
                                                                     "COLUMNCOUNT").getOutput(0)

    nrows = int(input_layer_prop.get("nrows", "!!!"))
    ncols = int(input_layer_prop.get("ncols", "!!!"))

    # try to get cell size set correctly for both x and y
    rowcs = float(input_layer_prop.get("cell size y", "!!!"))
    colcs = float(input_layer_prop.get("cell size x", "!!!"))

    def constrainsteps(n):
        total = 1
        while n >= 1500:
            n = (n // 2)
            total += 2
        total = total * 2
        return total, n

    nrow_steps, nrow_n = constrainsteps(nrows)
    ncol_steps, ncol_n = constrainsteps(ncols)

    cs_steps = ncol_steps * colcs
    rs_steps = nrow_steps * rowcs

    if not arcpy.Exists(os.path.join(output_directory, "Rasters", "tri")):
        arcpy.Resample_management(in_raster=input_tri,
                                  out_raster=os.path.join(output_directory, "Rasters", "tri"),
                                  cell_size=str(rs_steps) + " " + str(cs_steps), resampling_type="NEAREST")

    if not arcpy.Exists(os.path.join(output_directory, "Rasters", "bathy")):
        arcpy.Resample_management(in_raster=input_bathymetry,
                                  out_raster=os.path.join(output_directory, "Rasters", "bathy"),
                                  cell_size=str(rs_steps) + " " + str(cs_steps), resampling_type="NEAREST")

    inRas = arcpy.Raster(os.path.join(output_directory, "Rasters", "bathy"))
    xmin = inRas.extent.XMin
    ymin = inRas.extent.YMin
    xmax = inRas.extent.XMax
    ymax = inRas.extent.YMax

    nrows = int((inRas.extent.YMax - inRas.extent.YMin) / inRas.meanCellHeight)
    ncols = int((inRas.extent.XMax - inRas.extent.XMin) / inRas.meanCellWidth)

    del inRas

    fig = plt.figure()
    fig.subplots_adjust(hspace=0.1, wspace=0.1)

    arr = arcpy.da.FeatureClassToNumPyArray(os.path.join(output_directory, "Shapefiles", "CTD_location_points_p.shp"),
                                            ["SHAPE@X", "SHAPE@Y"])
    x, y = zip(*arr)

    for i in ["tri", "bathy"]:

        if i == "tri":
            grd = arcpy.RasterToNumPyArray(os.path.join(output_directory, "Rasters", i), arcpy.Point(xmin, ymin), ncols,
                                           nrows)
            grd = grd.astype(float)
            grd[grd == -9999.] = numpy.nan
            grd[grd == -99999.] = numpy.nan
            grd[grd < -3.4e+38] = numpy.nan
            ax = fig.add_subplot(1, 2, 1)
            ax.set_title("TRI (a)")
            im = ax.matshow(grd, extent=[xmin, xmax, ymin, ymax])
            plt.plot(x, y, 'ro', markersize=0)
        else:
            grd = arcpy.RasterToNumPyArray(os.path.join(output_directory, "Rasters", i), arcpy.Point(xmin, ymin), ncols,
                                           nrows)
            grd = grd.astype(float)
            grd[grd < -3.4e+38] = numpy.nan
            grd[grd == 3000000] = numpy.nan
            ax = fig.add_subplot(1, 2, 2)
            ax.set_title("Bathymetry (b)")
            im = ax.matshow(grd, extent=[xmin, xmax, ymin, ymax])
            plt.plot(x, y, 'ro', markersize=2)

        ax.set_xticks([])
        ax.set_yticks([])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        plt.colorbar(im, cax=cax)
        del grd, cax, divider
        gc.collect()
    del arr, ax
    gc.collect()

    fig.savefig(os.path.join(output_directory, "Plots", "tri_bathy.png"), dpi=300, bbox_inches='tight')
    fig.clf()
    plt.close()

    del fig, ax


def XYScatterPlot(x, y, plot_name):
    f, ax = plt.subplots()

    ax.plot(x, y, 'ro', markersize=4)
    # fig1 = lab.figure(figsize=(10.5, 10.5))
    # ax1 = fig1.add_axes([0.12, 0.15, 0.82, 0.75])

    ax.set_ylabel('Predicted (TRI/Bathymetry values)', multialignment='center', fontsize=16)
    ax.set_xlabel('Observed (CTD values)', multialignment='center', fontsize=16)

    # Get axis limits correct.
    if x.max() > y.max() and x.min() < y.min():
        ax.axis(ymin=x.min(), ymax=x.max(), xmin=x.min(), xmax=x.max())
    elif x.max() > y.max() and x.min() > y.min():
        ax.axis(ymin=y.min(), ymax=x.max(), xmin=y.min(), xmax=x.max())
    elif x.max() < y.max() and x.min() > y.min():
        ax.axis(ymin=y.min(), ymax=y.max(), xmin=y.min(), xmax=y.max())
    elif x.max() < y.max() and x.min() < y.min():
        ax.axis(ymin=x.min(), ymax=y.max(), xmin=x.min(), xmax=y.max())

    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")

    f.savefig(os.path.join(output_directory, "Plots", plot_name + ".png"), dpi=300, bbox_inches='tight')
    del f, ax


def BuildDFArcFC(table, index_col=None):
    # Useful function from: https://gist.github.com/om-henners/9b7008ac3e55b9e01eae97603a79c5f3
    """
    Build a pandas dataframe from an ArcGIS Feature Class.

    Uses the arcpy Search Cursor to loop over a feature class in pandas and
    generate a pandas dataframe. If the dataset is a feature class with a
    geometry it will calculate the length and the area before returning, and
    the geometry will be returned as well-known-text.

    :param table: The path to the feature class or table
    :type table: str
    :param index_col: A column to use as the dataframe index. If not supplied
                      for feature classes use the Object ID
    :type index_col: str
    :return: dataframe representation of the feature class. Note this is all in
             memory, so be careful with really big datasets!
    :rtype: pd.DataFrame
    """

    desc = arcpy.Describe(table)
    cursor = arcpy.SearchCursor(table)

    new_data = []

    for row in cursor:
        new_row = {}
        for field in desc.fields:
            new_row[field.aliasName or field.name] = row.getValue(field.name)
        new_data.append(new_row)

    try:
        if not index_col:
            index_col = desc.OIDFieldName

        df = pd.DataFrame(new_data).set_index(index_col)
        df["SHAPEArea"] = df[desc.shapeFieldName].apply(lambda g: g.area)
        df["SHAPELength"] = df[desc.shapeFieldName].apply(lambda g: g.length)
        df[desc.shapeFieldName] = df[desc.shapeFieldName].apply(lambda g: g.WKT)
    except AttributeError:
        # If this is a table in the datbase or on disk, in ArcGIS it won't have
        # either an OID field, nor a geometry
        pass

    return df


def CreateFishNet(input_bathymetry):
    arcpy.AddMessage("Now building fishnet... May take a short while..")

    if not arcpy.Exists(os.path.join(output_directory, "Rasters", "anull_rst")):
        arcpy.gp.SetNull_sa(input_bathymetry, input_bathymetry, os.path.join(output_directory, "Rasters", "anull_rst"),
                            "VALUE >= 0")
        input_bathymetry = os.path.join(output_directory, "Rasters", "anull_rst")

    if not arcpy.Exists(os.path.join(output_directory, "Rasters", "null_rst")):
        null_raster = Con(IsNull(input_bathymetry), 3000000, input_bathymetry)
        null_raster.save(os.path.join(output_directory, "Rasters", "null_rst"))

    input_bathymetry = os.path.join(output_directory, "Rasters", "null_rst")

    input_bathymetry_prop = {}

    input_bathymetry_prop["cell size x"] = arcpy.GetRasterProperties_management(input_bathymetry,
                                                                                "CELLSIZEX").getOutput(0)
    input_bathymetry_prop["cell size y"] = arcpy.GetRasterProperties_management(input_bathymetry,
                                                                                "CELLSIZEY").getOutput(0)
    input_bathymetry_prop["nrows"] = arcpy.GetRasterProperties_management(input_bathymetry,
                                                                          "ROWCOUNT").getOutput(0)
    input_bathymetry_prop["ncols"] = arcpy.GetRasterProperties_management(input_bathymetry,
                                                                          "COLUMNCOUNT").getOutput(0)
    input_bathymetry_desc = arcpy.Describe(input_bathymetry)
    input_bathymetry_prop["extent"] = input_bathymetry_desc.extent

    rExt = input_bathymetry_prop.get("extent", "!!!")

    # # try to get step size set correctly for both x and y
    # nrows = int(input_bathymetry_prop.get("nrows", "!!!"))
    # ncols = int(input_bathymetry_prop.get("ncols", "!!!"))
    #
    # # try to get cell size set correctly for both x and y
    # rowcs = float(input_bathymetry_prop.get("cell size y", "!!!"))
    # colcs = float(input_bathymetry_prop.get("cell size x", "!!!"))
    #
    # def constrainsteps(n):
    #     total = 1
    #     while n >= 200:
    #         n = (n // 2)
    #         total += 2
    #     total = total * 2
    #     return total, n
    #
    # nrow_steps, nrow_n = constrainsteps(nrows)
    # ncol_steps, ncol_n = constrainsteps(ncols)
    #
    # # arcpy.AddMessage("Resulting bathymetric tiles will be: " + str(nrow_n) + " by " + str(ncol_n) + " pixels.")
    #
    # xStep = rWid / ncol_steps  # the value of each step
    # yStep = rHgt / nrow_steps  # in the X and Y
    #
    # while xStep / colcs > 100.0:
    #     xStep = xStep / 2
    #     ncol_steps = ncol_steps * 2
    #
    # while yStep / rowcs > 100.0:
    #     yStep = yStep / 2
    #     nrow_steps = nrow_steps * 2
    #
    # list_of_sectors = []
    # for xIndex in range(ncol_steps):  # python iterating for var in list:
    #     for yIndex in range(nrow_steps):
    #         # calculate the extent, if you want overlap you can do it here
    #         Xmin = rExt.XMin + (xStep * xIndex)
    #         Xmax = Xmin + xStep
    #         Ymin = rExt.YMin + (yStep * yIndex)
    #         Ymax = Ymin + yStep
    #
    #         pol_array = arcpy.Array([arcpy.Point(Xmin, Ymin),
    #                                  arcpy.Point(Xmin, Ymax),
    #                                  arcpy.Point(Xmax, Ymax),
    #                                  arcpy.Point(Xmax, Ymin)])
    #
    #         rowcs_poly = Xmax - Xmin
    #         colcs_poly = Ymin - Ymax
    #
    #         list_of_sectors.append(arcpy.Polygon(pol_array))
    #
    # arcpy.CopyFeatures_management(list_of_sectors, os.path.join(output_directory, "Shapefiles", "fish_net_polygon.shp"))

    arcpy.CreateFishnet_management(
        out_feature_class=os.path.join(output_directory, "Shapefiles", "fish_net_polygon.shp"),
        origin_coord=str(rExt.XMin) + " " + str(rExt.YMin),
        y_axis_coord=str(rExt.XMin) + " " + str(rExt.YMax), cell_width="", cell_height="",
        number_rows="250", number_columns="250",
        corner_coord=str(rExt.XMax) + " " + str(rExt.YMax), labels="NO_LABELS",
        template=input_bathymetry,
        geometry_type="POLYGON")

    dsc = arcpy.Describe(input_bathymetry)
    coord_sys = dsc.spatialReference
    arcpy.DefineProjection_management(os.path.join(output_directory, "Shapefiles", "fish_net_polygon.shp"), coord_sys)

    arcpy.CalculateField_management(in_table=os.path.join(output_directory, "Shapefiles", "fish_net_polygon.shp"),
                                    field="Id", expression="[FID]",
                                    expression_type="VB", code_block="")

    rowcs_poly = (rExt.XMax - rExt.XMin) / 200

    arcpy.PolygonToRaster_conversion(in_features=os.path.join(output_directory, "Shapefiles", "fish_net_polygon.shp"),
                                     value_field="Id",
                                     out_rasterdataset=os.path.join(output_directory, "Rasters", "fishnet"),
                                     cell_assignment="CELL_CENTER", priority_field="NONE", cellsize=rowcs_poly)


def CalculateRMSE(input_predicted, input_observed):
    return ((input_predicted - input_observed) ** 2).mean() ** .5


def ExtractWOD(input_ctd, output_directory, input_bathymetry, input_tri):
    arcpy.AddMessage("Extracting bottom CTD values...")

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    if not os.path.exists(os.path.join(output_directory, "Shapefiles")):
        os.mkdir(os.path.join(output_directory, "Shapefiles"))
    if not os.path.exists(os.path.join(output_directory, "Rasters")):
        os.mkdir(os.path.join(output_directory, "Rasters"))
    if not os.path.exists(os.path.join(output_directory, "Plots")):
        os.mkdir(os.path.join(output_directory, "Plots"))

    if not os.path.exists(os.path.join(output_directory, "Shapefiles", 'CTD_locations_values_1.csv')):
        # build list of ctd files

        ctd_list = []

        for file in os.listdir(input_ctd):
            if file.endswith(".csv"):
                ctd_list.append(os.path.join(input_ctd, file))

        var_list = []

        for i in ctd_list:
            arcpy.AddMessage("... extracting variables from file: " + i)
            with open(i, 'rb') as csvfile:
                ctd_reader = csv.reader(csvfile, delimiter=',', quotechar='|')
                for row in ctd_reader:
                    if row[0].startswith("VARIABLES"):
                        for i in row:
                            if i.strip() not in var_list:
                                var_list.append(i.strip())

        # Clean up list
        var_list.remove("VARIABLES")
        var_list.remove("F")
        var_list.remove("O")
        var_list.remove("")
        var_list.append("cast")
        var_list.append("lat")
        var_list.append("lon")

        count_ctd_files = 1

        for i in ctd_list:
            arcpy.AddMessage("... extracting deepest CTD casts from file: " + i)

            with open(i, 'rb') as csvfile:
                ctd_reader = csv.reader(csvfile, delimiter=',', quotechar='|')

                need_header = True
                end_data = False
                run_once = True
                count = 0
                concat_df = []

                n_datasets = 0
                f_datasets = 0

                for row in ctd_reader:

                    if row[0].startswith("#--"):
                        if not need_header:
                            need_header = True
                        if end_data:
                            end_data = False

                    if row[0].startswith("END"):
                        count = 0
                        end_data = True

                    if need_header and not end_data:
                        for key in row:
                            if key.startswith("CAST"):
                                cast = row[2].strip(' ')
                            if key.startswith("Latitude"):
                                lat = row[2].strip(' ')
                            if key.startswith("Longitude"):
                                lon = row[2].strip(' ')
                            if key.startswith("VARIABLES"):
                                row_var = [x.strip(' ') for x in row]
                                row_var = row_var[:-1]
                            if key.startswith("Prof-Flag"):
                                need_header = False

                    elif not need_header and not end_data:
                        row = [x.strip(' ') for x in row]

                        if count == 0:
                            last_value = row[0]
                            previous_value = row[0]
                            count = count + 1
                        else:
                            last_value = row[0]

                        if last_value > previous_value:
                            last_row = row
                        else:
                            previous_value = row[0]

                    elif end_data:
                        try:
                            df = pd.DataFrame([last_row], columns=row_var)
                            df['cast'] = cast
                            df['lat'] = lat
                            df['lon'] = lon

                            row_var.remove("VARIABLES")
                            row_var = filter(lambda a: a != "F", row_var)
                            row_var = filter(lambda a: a != "O", row_var)
                            row_var = filter(lambda a: a != "", row_var)
                            row_var.append("cast")
                            row_var.append("lat")
                            row_var.append("lon")

                            df = df[row_var]

                            concat_df.append(df)
                            end_data = False

                            n_datasets = n_datasets + 1

                        except:
                            f_datasets = f_datasets + 1
                            pass

            concat_df_append = pd.concat(concat_df, ignore_index=True)
            concat_df_append.replace({'---0---': ''}, regex=True, inplace=True)

            concat_df_append = concat_df_append.sort('Depth', ascending=False).drop_duplicates(
                ['lat', 'lon']).sort_index()

            concat_df_append = concat_df_append.convert_objects(convert_numeric=True)

            concat_df_append['Depth'] = concat_df_append['Depth'] * -1

            concat_df_append.to_csv(
                os.path.join(output_directory, "Shapefiles", "CTD_locations_values_" + str(count_ctd_files) + ".csv"),
                index=False, header=True)
            count_ctd_files = count_ctd_files + 1
            del df, row_var
    else:
        ctd_list = []
        for file in os.listdir(os.path.join(output_directory, "Shapefiles")):
            if file.endswith(".csv"):
                ctd_list.append(os.path.join(output_directory, "Shapefiles", file))

        frame = pd.DataFrame()
        list_ = []
        for file_ in ctd_list:
            df = pd.read_csv(file_, index_col=None, header=0)
            list_.append(df)

        concat_df_append = pd.concat(list_)

    # CREATE SHAPEFILE OF POINTS WITH COORDS OF CTD CASTS.
    if not os.path.exists(os.path.join(output_directory, "Shapefiles", "CTD_location_points.shp")):

        output_ctd_shp = os.path.join(output_directory, "Shapefiles", "CTD_location_points.shp")

        concat_df_append = concat_df_append.fillna(-999999)

        location_array = numpy.array(concat_df_append.to_records())

        arcpy.da.NumPyArrayToFeatureClass(location_array, output_ctd_shp, ['lon', 'lat'], arcpy.SpatialReference(4326))

        if not arcpy.Exists(os.path.join(output_directory, "Shapefiles", "CTD_location_points_p.shp")):
            arcpy.Project_management(in_dataset=os.path.join(output_directory, "Shapefiles", "CTD_location_points.shp"),
                                     out_dataset=os.path.join(output_directory, "Shapefiles",
                                                              "CTD_location_points_p.shp"),
                                     out_coor_system=input_bathymetry,
                                     transform_method="",
                                     in_coor_system=os.path.join(output_directory, "Shapefiles",
                                                                 "CTD_location_points.shp"),
                                     preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")

        output_ctd_shp = os.path.join(output_directory, "Shapefiles", "CTD_location_points_p.shp")

        ExtractMultiValuesToPoints(output_ctd_shp, [[input_bathymetry, "bathy_dep"]], "NONE")

        with arcpy.da.UpdateCursor(output_ctd_shp, "bathy_dep") as cursor:
            for row in cursor:
                if row[0] == -9999:
                    cursor.deleteRow()

        ExtractMultiValuesToPoints(output_ctd_shp, [[input_tri, "tri"]], "NONE")

        del location_array, concat_df_append

    if not os.path.exists(os.path.join(output_directory, "Shapefiles", "fish_net_polygon.shp")):
        CreateFishNet(input_bathymetry)

    if not os.path.exists(os.path.join(output_directory, "Plots", "tri_bathy.png")):
        ArcpyMap_1(os.path.join(output_directory, "Rasters", "null_rst"), input_tri, output_directory)

    if not os.path.exists(os.path.join(output_directory, "Shapefiles", "CTD_location_points_p_sj.shp")):
        arcpy.SpatialJoin_analysis(
            target_features=os.path.join(output_directory, "Shapefiles", "CTD_location_points_p.shp"),
            join_features=os.path.join(output_directory, "Shapefiles", "fish_net_polygon.shp"),
            out_feature_class=os.path.join(output_directory, "Shapefiles", "CTD_location_points_p_sj.shp"),
            join_operation="JOIN_ONE_TO_MANY", join_type="KEEP_ALL", field_mapping="",
            match_option="WITHIN", search_radius="", distance_field_name="")


def ValidateLayerStatisticsPlots(output_directory, tri_variable, depth_cut_off):
    arcpy.AddMessage("Calculating statistics and generating plots...")

    if tri_variable == "temp":
        wod_var_name = "Temperatur"
        woa_units = "degrees C"
    elif tri_variable == "sal":
        wod_var_name = "Salinity"
        woa_units = "unitless"
    elif tri_variable == "diso2":
        wod_var_name = "Oxygen"
        woa_units = "ml per l"
    else:
        arcpy.AddMessage("No matching TRI variable name, unable to proceed, check value for 'tri_variable'..")
        sys.exit(0)

    location_array = BuildDFArcFC(os.path.join(output_directory, "Shapefiles", "CTD_location_points_p_sj.shp"))
    location_array = location_array[location_array[wod_var_name] != -999999]
    location_array = location_array[location_array[wod_var_name] != -9999]

    # Now we need to test the bathymetry against the CTD data (we will say if <> x% difference, we scrap point)
    location_array['dep_diff'] = location_array['Depth'] / location_array['bathy_dep']

    location_array_0 = location_array[location_array['dep_diff'] >= 1.0 - depth_cut_off]
    location_array_1 = location_array_0[location_array_0['dep_diff'] <= 1.0 + depth_cut_off]

    n_dropped_ctd = len(location_array.index) - len(location_array_1.index)
    n_retained = len(location_array_1.index)

    rms_original_dep = CalculateRMSE(location_array['bathy_dep'], location_array['Depth'])
    rms_filtered_dep = CalculateRMSE(location_array_1['bathy_dep'], location_array_1['Depth'])

    rms_original_tri = CalculateRMSE(location_array['tri'], location_array[wod_var_name])
    rms_filtered_tri = CalculateRMSE(location_array_1['tri'], location_array_1[wod_var_name])

    location_array['tri_diff'] = location_array[wod_var_name] - location_array['tri']
    location_array_1['tri_diff'] = location_array_1[wod_var_name] - location_array_1['tri']

    diff_original_tri = location_array['tri_diff'].mean()
    diff_filtered_tri = location_array_1['tri_diff'].mean()

    corr_pearson_original_tri, corr_pearson_p_original_tri = stats.pearsonr(location_array[wod_var_name],
                                                                            location_array['tri'])
    corr_pearson_filtered_tri, corr_pearson_p_filtered_tri = stats.pearsonr(location_array_1[wod_var_name],
                                                                            location_array_1['tri'])

    # Figure Depths - XY Scatter plots
    if not os.path.exists(os.path.join(output_directory, "Plots", "original_depth.png")):
        XYScatterPlot(location_array['Depth'], location_array['bathy_dep'], "original_depth")
        XYScatterPlot(location_array_1['Depth'], location_array_1['bathy_dep'], "filtered_depth")
        XYScatterPlot(location_array[wod_var_name], location_array['tri'], "original_tri")
        XYScatterPlot(location_array_1[wod_var_name], location_array_1['tri'], "filtered_tri")

    # Depth bin error plot
    if not os.path.exists(os.path.join(output_directory, "Plots", "tri_by_depth.png")):
        bins = [-15000, -10000, -8000, -6000, -5000, -4000, -3000, -2000, -1500, -1250, -1000, -750, -500, -300, -200,
                -100, -50, 0]
        labels = [-10000, -8000, -6000, -5000, -4000, -3000, -2000, -1500, -1250, -1000, -750, -500, -300, -200, -100,
                  -50, 0]

        lists = [[]]
        names = []
        counts = []

        location_array_1['binned'] = pd.cut(location_array_1['Depth'], bins=bins, labels=labels)

        for name, group in location_array_1.groupby('binned'):
            names.append(name)
            lists.append(group['tri_diff'])
            counts.append(len(group))
        f, ax = plt.subplots()

        counts = [0] + counts
        names = [-15000] + names
        counts_names = zip(names, counts)

        bp = ax.boxplot(lists, patch_artist=True)
        ax.set_xticklabels(['%s (%d)' % (k, (v)) for k, v in counts_names], rotation='vertical')
        ax.set_ylabel('Observed cast value - Predicted TRI', multialignment='center', fontsize=16)
        ax.set_xlabel('Depth (m)', multialignment='center', fontsize=16)
        ax.axhline(y=0, color='gray', linestyle='--')

        for box in bp['boxes']:
            box.set(color='#7570b3', linewidth=1)
            box.set(facecolor='#50A6C2')

        f.savefig(os.path.join(output_directory, "Plots", "tri_by_depth.png"), dpi=300, bbox_inches='tight')
        del f, ax, bp, names, lists, name, group

    # Spatial error plot
    SpatialErrorMap(depth_cut_off, wod_var_name)

    # Generate HTML output file with various details
    html_str = "<html><head><title>GlobENV - %s</title></head><body><h1>GlobENV Validation Output - %s</h1>" \
               "<p>This page contains some analysis for the bathymetry and generated " \
               "trilinearly interpolated variable: <b>%s</b>. If you want to conduct further " \
               "analyses, a list of generated files can be found at the end of this file.</p><hr>" % (
               tri_variable, tri_variable, tri_variable)

    html_str += "<h2>Overall pattern</h2>The following map shows: The output TRI variable (a), the depth layer " \
                "with locations of the CTD points as red points (b). You should look at where these points fall, " \
                "large clusters, especially in shallower water areas or large expanses of relatively flat " \
                "bathymetry  will generally explain poor model performance." \
                "<p><img src=\"%s\" width=\"100%%\"></p><hr>" % os.path.join("Plots", "tri_bathy.png")

    html_str += "<h2>TRI Model Performance</h2>The following plots and statistics test the performance of the " \
                "model against validation data that you downloaded from World Ocean database. <p><b>Unfiltered data:</b> In the plots below, the first column shows " \
                "the comparison between the estimated bottom depth from the deepest CTD cast, and the second, the comparison of " \
                "the TRI model prior to filtering the data for casts that were less than %d %% of the bathymetric value for " \
                "the location of the cast.</p>" \
                "<p><table width=\"100%%\"><tr><td><img src=\"%s\" width=\"100%%\"></td><td><img src=\"%s\" width=\"100%%\"></td></tr>" \
                "</table></p>" % (depth_cut_off * 100, os.path.join("Plots", "original_depth.png"),
                                  os.path.join("Plots", "original_tri.png"))

    html_str += "<p><b>Filtered data:</b> The first column shows the comparison between the estimated bottom depth from the deepest CTD cast, " \
                "and the second, the comparison of the TRI model after filtering the data for casts that were within %d %% of the bathymetric value for " \
                "the location of the cast. " \
                "<p><table width=\"100%%\"><tr><td><img src=\"%s\" width=\"100%%\"></td><td><img src=\"%s\" width=\"100%%\"></td></tr>" \
                "</table></p>" % (depth_cut_off * 100, os.path.join("Plots", "filtered_depth.png"),
                                  os.path.join("Plots", "filtered_tri.png"))

    html_str += "<p><b>Performance at different depths:</b> This plot shows how the TRI variable performs at different depths. Often, the TRI approach " \
                "does not perform well in shallow waters or large geographical areas that change little with depth and transition across multiple watermasses. " \
                "Ideally, the best model will have low difference across the whole of the depth range. If the values are above the zero line, then the validation " \
                "CTD casts are less than the predicted TRI for the location, below it, TRI is lower." \
                "<p><img src=\"%s\" width=\"50%%\"></p>" % os.path.join("Plots", "tri_by_depth.png")

    html_str += "<p><b>Performance at different locations:</b> These maps show how the TRI variable performs at different locations. The first map (a) shows the " \
                "number of casts within a given validation cell. The second shows the difference between the observed CTD data and the predicted TRI model. If the " \
                "values are above the zero line, then the validation CTD casts are less than the predicted TRI for the location, below it, TRI is lower." \
                "<p><img src=\"%s\" width=\"100%%\"></p>" % os.path.join("Plots", "spatial_error.png")

    html_str += "<hr><h2>Summary statistics</h2><p>The statistics provide an insight into the overall performance of the " \
                "model against validation data that you downloaded from World Ocean database.</p>"

    html_str += """
                <table border="1">
                  <tr><b>
                    <th>Test</th>
                    <th>Unfiltered</th>
                    <th>Filtered</th>
                    <th>Notes</th>
                  </tr></b>
                   <tr>
                    <td>Samples</td>
                    <td>%s</td>
                    <td>%s</td>
                    <td>The number of CTD casts used for validation.</td>
                  </tr>
                  <tr>
                    <td>RMSE (Depth)</td>
                    <td>%s</td>
                    <td>%s</td>
                    <td>Root mean square error, this is overall error of the depth layer.</td>
                  </tr>
                  <tr>
                    <td>RMSE (TRI)</td>
                    <td>%s</td>
                    <td>%s</td>
                    <td>Root mean square error, this is overall error of the tri layer.</td>
                  </tr>
                  <tr>
                    <td>Mean difference (Cast - TRI)</td>
                    <td>%s</td>
                    <td>%s</td>
                    <td>Mean overall difference between validation cast and the modelled TRI layer.</td>
                  </tr>
                  <td>Pearson Correlation (Cast vs TRI)</td>
                    <td>%s</td>
                    <td>%s</td>
                    <td>Pearson correlation of the validation casts against modelled TRI per location.</td>
                  </tr>
                </table>
                """ % (n_dropped_ctd + n_retained, n_retained, rms_original_dep, rms_filtered_dep, rms_original_tri,
                       rms_filtered_tri, diff_original_tri, diff_filtered_tri,
                       str(round(corr_pearson_original_tri, 2)) + " (p = " +
                       str(round(corr_pearson_p_original_tri, 2)) + ")",
                       str(round(corr_pearson_filtered_tri, 2)) + " (p = " +
                       str(round(corr_pearson_p_filtered_tri, 2)) + ")")

    html_str += "<hr><h2>Output files</h2><p>Generated output files are as follows:</p> " \
                "<ol><li>Spatially joined CTD location points: %s </li>" \
                "<li>Fishnet grid: %s </li>" \
                "<li>All Plots are in: %s </li></ol>" % (
                os.path.join(output_directory, "Shapefiles", "CTD_location_points_p_sj.shp"),
                os.path.join(output_directory, "Shapefiles", "fish_net_polygon.shp"),
                os.path.join(output_directory, "Plots"))

    html_file = open(os.path.join(output_directory, "validation.html"), "w")
    html_file.write(html_str)
    html_file.close()


if __name__ == "__main__":
    # input ctd data directory (these are the downloaded unzipped files that come from WOD).
    input_ctd = r"E:\2_Bathymetries\Validation"
    # set an output directory, best to be empty
    output_directory = r"E:\2_Bathymetries\Validation\Output_Sal"
    # your original bathymetric file that you used for the bathymetry data preparation script
    input_bathymetry = r"E:\2_Bathymetries\Mediterranean\Layers\depth.tif"
    # input tri
    input_tri = r"E:\2_Bathymetries\Mediterranean\Layers\sal.tif"
    # variable to process ("temp", "salinity")
    tri_variable = "sal"
    # depth cut off - if CTD max depth is proportionally different than bathymetric layer depth
    depth_cut_off = 0.1  # 0.1 = 10% plus and minus

    ExtractWOD(input_ctd, output_directory, input_bathymetry, input_tri)
    ValidateLayerStatisticsPlots(output_directory, tri_variable, depth_cut_off)

    # # input ctd data directory (these are the downloaded unzipped files that come from WOD).
    # input_ctd = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\ETOPO2\Validation"
    # # set an output directory, best to be empty
    # output_directory = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco08\Validation\sal"
    # # your original bathymetric file that you used for the bathymetry data preparation script
    # input_bathymetry = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco08\Bathymetry\gebco08p"
    # # input tri
    # input_tri = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco08\Outputs\sal"
    # # variable to process ("temp", "salinity")
    # tri_variable = "sal"
    # # depth cut off - if CTD max depth is proportionally different than bathymetric layer depth
    # depth_cut_off = 0.1  # 0.1 = 10% plus and minus
    #
    # ExtractWOD(input_ctd, output_directory, input_bathymetry, input_tri)
    # ValidateLayerStatisticsPlots(output_directory, tri_variable, depth_cut_off)

    # # input ctd data directory (these are the downloaded unzipped files that come from WOD).
    # input_ctd = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\ETOPO2\Validation"
    # # set an output directory, best to be empty
    # output_directory = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco08\Validation\diso2"
    # # your original bathymetric file that you used for the bathymetry data preparation script
    # input_bathymetry = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco08\Bathymetry\gebco08p"
    # # input tri
    # input_tri = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco08\Outputs\diso2"
    # # variable to process ("temp", "salinity")
    # tri_variable = "diso2"
    # # depth cut off - if CTD max depth is proportionally different than bathymetric layer depth
    # depth_cut_off = 0.1  # 0.1 = 10% plus and minus
    #
    # ExtractWOD(input_ctd, output_directory, input_bathymetry, input_tri)
    # ValidateLayerStatisticsPlots(output_directory, tri_variable, depth_cut_off)

    # gc.collect()
    #
    # # input ctd data directory (these are the downloaded unzipped files that come from WOD).
    # input_ctd = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\ETOPO2\Validation"
    # # set an output directory, best to be empty
    # output_directory = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco14\Validation\temp"
    # # your original bathymetric file that you used for the bathymetry data preparation script
    # input_bathymetry = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco14\Bathymetry\gebco14p"
    # # input tri
    # input_tri = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco14\Outputs\temp"
    # # variable to process ("temp", "salinity")
    # tri_variable = "temp"
    # # depth cut off - if CTD max depth is proportionally different than bathymetric layer depth
    # depth_cut_off = 0.1  # 0.1 = 10% plus and minus
    #
    # ExtractWOD(input_ctd, output_directory, input_bathymetry, input_tri)
    # ValidateLayerStatisticsPlots(output_directory, tri_variable, depth_cut_off)
    # gc.collect()
    # # input ctd data directory (these are the downloaded unzipped files that come from WOD).
    # input_ctd = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\ETOPO2\Validation"
    # # set an output directory, best to be empty
    # output_directory = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco14\Validation\sal"
    # # your original bathymetric file that you used for the bathymetry data preparation script
    # input_bathymetry = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco14\Bathymetry\gebco14p"
    # # input tri
    # input_tri = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco14\Outputs\sal"
    # # variable to process ("temp", "salinity")
    # tri_variable = "sal"
    # # depth cut off - if CTD max depth is proportionally different than bathymetric layer depth
    # depth_cut_off = 0.1  # 0.1 = 10% plus and minus
    #
    # ExtractWOD(input_ctd, output_directory, input_bathymetry, input_tri)
    # ValidateLayerStatisticsPlots(output_directory, tri_variable, depth_cut_off)
    # gc.collect()
    # # input ctd data directory (these are the downloaded unzipped files that come from WOD).
    # input_ctd = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\ETOPO2\Validation"
    # # set an output directory, best to be empty
    # output_directory = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco14\Validation\diso2"
    # # your original bathymetric file that you used for the bathymetry data preparation script
    # input_bathymetry = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco14\Bathymetry\gebco14p"
    # # input tri
    # input_tri = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco14\Outputs\diso2"
    # # variable to process ("temp", "salinity")
    # tri_variable = "diso2"
    # # depth cut off - if CTD max depth is proportionally different than bathymetric layer depth
    # depth_cut_off = 0.1  # 0.1 = 10% plus and minus
    #
    # ExtractWOD(input_ctd, output_directory, input_bathymetry, input_tri)
    # ValidateLayerStatisticsPlots(output_directory, tri_variable, depth_cut_off)
    # gc.collect()
    #
    # # input ctd data directory (these are the downloaded unzipped files that come from WOD).
    # input_ctd = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\ETOPO2\Validation"
    # # set an output directory, best to be empty
    # output_directory = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\srtm30\Validation\temp"
    # # your original bathymetric file that you used for the bathymetry data preparation script
    # input_bathymetry = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\srtm30\Bathymetry\srtm30p"
    # # input tri
    # input_tri = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\srtm30\Outputs\temp"
    # # variable to process ("temp", "salinity")
    # tri_variable = "temp"
    # # depth cut off - if CTD max depth is proportionally different than bathymetric layer depth
    # depth_cut_off = 0.1  # 0.1 = 10% plus and minus
    #
    # ExtractWOD(input_ctd, output_directory, input_bathymetry, input_tri)
    # ValidateLayerStatisticsPlots(output_directory, tri_variable, depth_cut_off)
    # gc.collect()
    # # input ctd data directory (these are the downloaded unzipped files that come from WOD).
    # input_ctd = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\ETOPO2\Validation"
    # # set an output directory, best to be empty
    # output_directory = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\srtm30\Validation\sal"
    # # your original bathymetric file that you used for the bathymetry data preparation script
    # input_bathymetry = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\srtm30\Bathymetry\srtm30p"
    # # input tri
    # input_tri = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\srtm30\Outputs\sal"
    # # variable to process ("temp", "salinity")
    # tri_variable = "sal"
    # # depth cut off - if CTD max depth is proportionally different than bathymetric layer depth
    # depth_cut_off = 0.1  # 0.1 = 10% plus and minus
    #
    # ExtractWOD(input_ctd, output_directory, input_bathymetry, input_tri)
    # ValidateLayerStatisticsPlots(output_directory, tri_variable, depth_cut_off)
    # gc.collect()
    # # input ctd data directory (these are the downloaded unzipped files that come from WOD).
    # input_ctd = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\ETOPO2\Validation"
    # # set an output directory, best to be empty
    # output_directory = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\srtm30\Validation\diso2"
    # # your original bathymetric file that you used for the bathymetry data preparation script
    # input_bathymetry = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\srtm30\Bathymetry\srtm30p"
    # # input tri
    # input_tri = r"//10.10.10.104\Model_Outputs\GlobEnv\2_Bathymetries_and_Outputs\Global\srtm30\Outputs\diso2"
    # # variable to process ("temp", "salinity")
    # tri_variable = "diso2"
    # # depth cut off - if CTD max depth is proportionally different than bathymetric layer depth
    # depth_cut_off = 0.1  # 0.1 = 10% plus and minus
    #
    # ExtractWOD(input_ctd, output_directory, input_bathymetry, input_tri)
    # ValidateLayerStatisticsPlots(output_directory, tri_variable, depth_cut_off)
    # gc.collect()