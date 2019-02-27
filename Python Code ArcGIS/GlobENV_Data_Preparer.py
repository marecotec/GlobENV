#  GlobENV Data Preparation

import arcpy
import os, csv
from arcpy import env
from Includes import load_depth_string

import pandas as pd
import gc
import numpy as np
import time
import multiprocessing
from functools import partial
import shutil
import sys

gc.enable()

arcpy.CheckOutExtension("Spatial")


# This function pre-reads data to assess whether or not this script
# can handle it, and prepare it for GlobENV

def data_check(output_directory, input_netcdf):

    # Steps
    # 1 - Try to read NetCDF
    # 2 - Is NetCDF in appropriate 3D format?
    # 3 - What is the vertical grid separation?
    # 4 - Which variables are in the dataset?
    # 5 - Select the variables etc that you need interactively.
    # 6 - Print out a config file that can be used

    try:
        # 1 - Try to read NetCDF
        if input_netcdf is not None:
            arcpy.AddMessage("Attempting to reading file: " + str(input_netcdf))
            if NetCDFFile.isNetCDF(input_netcdf):
                arcpy.AddMessage(".. file read successful")
                netCDFFile = NetCDFFile(input_netcdf)
                netcdf_variables = netCDFFile.getVariables()
            else:
                arcpy.AddMessage(".. file read unsuccessful, check file exists or interrogating file in Panoply.")

        





    except:
        print "test"




def mpprocess(output_directory, variable_name, input_woa_netcdf,
              interpolation_procedure, interpolation_resolution,
              coordinate_system, extraction_extent,
              createxyz, depth_range):

    try:
        i = depth_range

        status = "Started"

        arcpy.env.extent = extraction_extent

        arcpy.AddMessage("Processing depth: " + str(int(i)))

        out_temp_layer = variable_name[0:4] + str(int(i)) + ".shp"
        dimensionValues = "depth " + str(int(i))

        output_dir_geographic = os.path.join(output_directory, "temp", "Geographic",
                                             variable_name[0:4] + str(int(i)))
        output_geographic = os.path.join(output_dir_geographic, variable_name[0:4] + str(int(i)))

        output_dir_projected = os.path.join(output_directory, "temp", "Projected", variable_name[0:4] + str(int(i)))
        output_projected = os.path.join(output_dir_projected, variable_name[0:4] + str(int(i)))

        if not os.path.exists(os.path.join(output_dir_projected, "xy_coords.yxz")) or os.path.exists(
                os.path.join(output_dir_geographic, "xy_coords.yxz")):

            if not os.path.exists(output_dir_geographic):
                os.makedirs(output_dir_geographic)

            if not os.path.exists(output_dir_projected):
                os.makedirs(output_dir_projected)

            env.workspace = os.path.join(output_directory, "temp", variable_name[0:4] + str(int(i)))

            # 1 Extract layer to a temporary feature class
            arcpy.MakeNetCDFFeatureLayer_md(in_netCDF_file=input_woa_netcdf, variable=variable_name,
                                            x_variable="lon",
                                            y_variable="lat",
                                            out_feature_layer=out_temp_layer,
                                            row_dimension="lat;lon",
                                            z_variable="", m_variable="", dimension_values=dimensionValues,
                                            value_selection_method="BY_VALUE")

            # 2 Interpolate to higher resolution and 3 save to output directory
            if interpolation_procedure == "IDW":
                status = "Interpolating " + str(int(i)) + " using IDW"
                arcpy.gp.Idw_sa(out_temp_layer, variable_name, output_geographic,
                                interpolation_resolution, "2", "VARIABLE 10", "")
            elif interpolation_procedure == "Spline":
                status = "Interpolating " + str(int(i)) + " using Spline"
                arcpy.CopyFeatures_management(out_temp_layer, os.path.join(output_dir_geographic, "out.shp"))
                arcpy.gp.Spline_sa(out_temp_layer, variable_name, output_geographic,
                                   interpolation_resolution, "TENSION", "0.1", "10")
                arcpy.Delete_management(
                    os.path.join(output_directory, "temp", "Geographic", variable_name[0:4] + str(int(i)),
                                 "out.shp"))
            elif interpolation_procedure == "Kriging":
                status = "Interpolating " + str(int(i)) + " using Ordinary Kriging"
                arcpy.gp.Kriging_sa(out_temp_layer, variable_name,
                                    output_geographic,
                                    "Spherical " + str(interpolation_resolution), interpolation_resolution,
                                    "VARIABLE 10", "")

            elif interpolation_procedure == "None":
                status = "Making a raster for " + str(int(i))
                arcpy.MakeNetCDFRasterLayer_md(in_netCDF_file=input_woa_netcdf, variable=variable_name,
                                               x_dimension="lon", y_dimension="lat",
                                               out_raster_layer=output_geographic,
                                               band_dimension="", dimension_values="",
                                               value_selection_method="BY_VALUE")

            if len(coordinate_system) > 1:
                status = "Reprojecting " + variable_name[0:4] + str(int(i)) + "."
                arcpy.ProjectRaster_management(output_geographic,
                                               output_projected,
                                               coordinate_system, "NEAREST", "#", "#", "#", "#")

            if createxyz == "Only Geographic" or createxyz == "Both":
                status = "Building geographic xy coords for " + variable_name[0:4] + str(int(i)) + "."
                raster_to_xyz(output_geographic,
                              variable_name[0:4] + str(int(i)),
                              output_dir_geographic,
                              349000000.0)

                df = pd.read_csv(os.path.join(output_dir_geographic, variable_name[0:4] + str(int(i)) + ".yxz"),
                                 header=0, names=["y", "x", "z"], sep=" ",
                                 dtype={"y": np.float32, "x": np.float32, "z": np.float32})
                master = df[["x", "y"]].copy()
                master.columns = ["x", "y"]
                master = np.round(master, 4)
                master.to_pickle(os.path.join(output_dir_geographic, "xy_coords.yxz"))
                del master, df
                gc.collect()

            if createxyz == "Only Projected" or createxyz == "Both":
                status = "Building projected xy coords for " + variable_name[0:4] + str(int(i)) + "."
                raster_to_xyz(output_projected,
                              variable_name[0:4] + str(int(i)),
                              output_dir_projected,
                              349000000.0)

                df = pd.read_csv(os.path.join(output_dir_projected, variable_name[0:4] + str(int(i)) + ".yxz"),
                                 header=0, names=["y", "x", "z"], sep=" ",
                                 dtype={"y": np.float32, "x": np.float32, "z": np.float32})
                master = df[["x", "y"]].copy()
                master.columns = ["x", "y"]
                master = np.round(master, 4)
                master.to_pickle(os.path.join(output_dir_projected, "xy_coords.yxz"))
                del master, df
                gc.collect()
        else:
            arcpy.AddMessage("Skipping " + str(int(i)) + ".")
    except:
        arcpy.AddMessage("Failed on " + str(int(i)) + ",at status: " + str(status))
        arcpy.AddMessage(arcpy.GetMessages())

    def execute(self, parameters, messages):

        t_start = time.clock()

        arcpy.env.overwriteOutput = True

        for param in parameters:
            arcpy.AddMessage("Parameter: %s = %s" % (param.name, param.valueAsText))

        input_woa_netcdf = parameters[0].valueAsText
        variable_name = parameters[1].valueAsText
        depths = parameters[2].valueAsText
        interpolation_procedure = parameters[3].valueAsText
        interpolation_resolution = parameters[4].valueAsText
        extraction_extent = parameters[5].valueAsText
        output_directory = parameters[6].valueAsText
        coordinate_system = parameters[7].valueAsText
        createxyz = parameters[8].valueAsText
        cpu_cores_used = parameters[9].valueAsText

        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        if not os.path.exists(os.path.join(output_directory, "Projected")):
            os.makedirs(os.path.join(output_directory, "Projected"))

        if not os.path.exists(os.path.join(output_directory, "Geographic")):
            os.makedirs(os.path.join(output_directory, "Geographic"))

        arcpy.env.extent = extraction_extent

        arcpy.AddMessage("Extracting " + str(input_woa_netcdf) + ".")

        # Set environment variables and build other variables for processing
        arcpy.env.mask = ""
        arcpy.env.workspace = output_directory
        depth_range = load_depth_string(depths)

        arcpy.AddMessage("There are " + str(len(depth_range)) + " to process.")

        if int(cpu_cores_used) > int(multiprocessing.cpu_count()):
            cpu_cores_used = multiprocessing.cpu_count() - 1

        arcpy.AddMessage("Will use " + str(cpu_cores_used) + " cores for processing")

        python_exe = os.path.join(sys.exec_prefix, 'pythonw.exe')
        multiprocessing.set_executable(python_exe)

        pool = multiprocessing.Pool(int(cpu_cores_used))
        func = partial(mpprocess_call, output_directory, variable_name, input_woa_netcdf, interpolation_procedure,
                       interpolation_resolution, coordinate_system, extraction_extent, createxyz)
        pool.map(func, depth_range)
        pool.close()
        pool.join()

        count = 0

        depth_range = load_depth_string(depths)

        for i in depth_range:
            try:
                output_geographic = os.path.join(output_directory, "temp", "Geographic",
                                                 variable_name[0:4] + str(int(i)), variable_name[0:4] + str(int(i)))
                if arcpy.Exists(output_geographic) and not arcpy.Exists(
                        os.path.join(output_directory, "Geographic", variable_name[0:4] + str(int(i)))):
                    shutil.copytree(output_geographic, os.path.join(output_directory, "Geographic",
                                                                    variable_name[0:4].lower() + str(int(i))))
            except:
                arcpy.AddMessage("Issue copying, geographic for depth " + str(int(i)))

            try:
                output_projected = os.path.join(output_directory, "temp", "Projected", variable_name[0:4] + str(int(i)),
                                                variable_name[0:4] + str(int(i)))
                if arcpy.Exists(output_projected) and not arcpy.Exists(
                        os.path.join(output_directory, "Projected", variable_name[0:4] + str(int(i)))):
                    shutil.copytree(output_projected, os.path.join(output_directory, "Projected",
                                                                   str(variable_name[0:4].lower() + str(int(i)))))
            except:
                arcpy.AddMessage("Issue copying, projected for depth " + str(int(i)))

            try:
                if count == 0:
                    if os.path.exists(
                            os.path.join(output_directory, "temp", "Projected", variable_name[0:4] + str(int(i)),
                                         "xy_coords.yxz")):
                        shutil.copyfile(
                            os.path.join(output_directory, "temp", "Projected", variable_name[0:4] + str(int(i)),
                                         "xy_coords.yxz"),
                            os.path.join(output_directory, "Projected", "xy_coords.yxz"))
                        count = 1
            except:
                arcpy.AddMessage("Issue copying, xyz coordiantes for depth " + str(int(i)))

        arcpy.AddMessage("Making pyramids and statistics for outputs")
        arcpy.BuildPyramidsandStatistics_management(in_workspace=os.path.join(output_directory, "Geographic"),
                                                    include_subdirectories="NONE",
                                                    build_pyramids="BUILD_PYRAMIDS",
                                                    calculate_statistics="CALCULATE_STATISTICS", BUILD_ON_SOURCE="NONE",
                                                    block_field="", estimate_statistics="NONE", x_skip_factor="1",
                                                    y_skip_factor="1", ignore_values="", pyramid_level="-1",
                                                    SKIP_FIRST="NONE", resample_technique="NEAREST",
                                                    compression_type="DEFAULT", compression_quality="75",
                                                    skip_existing="SKIP_EXISTING")
        if len(coordinate_system) > 1:
            arcpy.BuildPyramidsandStatistics_management(in_workspace=os.path.join(output_directory, "Projected"),
                                                        include_subdirectories="NONE",
                                                        build_pyramids="BUILD_PYRAMIDS",
                                                        calculate_statistics="CALCULATE_STATISTICS",
                                                        BUILD_ON_SOURCE="NONE",
                                                        block_field="", estimate_statistics="NONE", x_skip_factor="1",
                                                        y_skip_factor="1", ignore_values="", pyramid_level="-1",
                                                        SKIP_FIRST="NONE", resample_technique="NEAREST",
                                                        compression_type="DEFAULT", compression_quality="75",
                                                        skip_existing="SKIP_EXISTING")

        arcpy.AddMessage("Script complete in %s hours." % str((time.clock() - t_start) / 3600))

        return


class NetCDFFile(object):
    """-------------------------------------------------------------------------------------
    Class Name: NetCDFFile
    Creation Date: 3/1/2012
    Creator: KSigwart
    Description: This is a helper class that extends the NetCDFFile Properties to be more in
                 line with tools that let users interact with a netCDF file in forms of
                 lat/lon points.
    Inputs:
            netCDF file path:  The location of a netCDF file
    -------------------------------------------------------------------------------------"""

    # Class Parameters
    __sourceLocation = ''
    __variables = list()
    __dimensions = list()
    __latDimension = 'lat'
    __lonDimension = 'lon'
    __latValue = 0
    __lonValue = 0

    def __init__(self, netCDFloc):
        '''
         Defines the Class Properties based off the the NetCDF File Properties inputted
        '''
        ncFileProp = arcpy.NetCDFFileProperties(netCDFloc)
        self.__ncFileProperties = ncFileProp

        self.__sourceLocation = netCDFloc;
        self.__determineDimensions()
        self.__determineVariables()

    def __determineDimensions(self):
        '''
         We only want the dimensions that are not lat and lon.  The lat, lon
         dimensions are already being described by the point that will be mapped
         by the user.
        '''
        dimensions = self.__ncFileProperties.getDimensions();

        # Storing lat and lon dimensions to get the varables associated with both
        dimList = list(dimensions)
        newDimList = list()
        for dim in dimList:
            # if str(dim).lower().startswith('lat'):
            if 'lat' in str(dim).lower():
                self.__latDimension = str(dim)
            # elif str(dim).lower().startswith('lon'):
            elif 'lon' in str(dim).lower():
                self.__lonDimension = str(dim)
            else:
                newDimList.append(dim);

        self.__dimensions = newDimList
        return newDimList

    def __determineVariables(self):
        '''
         We only want variables that have a lat and long dimensions b/c we are
         mapping the variable to a point.  Also, lat and lon are already described
         by the point.
        '''
        latVariables = list(self.__ncFileProperties.getVariablesByDimension(self.__latDimension))
        lonVariables = list(self.__ncFileProperties.getVariablesByDimension(self.__lonDimension))

        variables = list()

        for variable in latVariables:
            if (variable != self.__latDimension and variable != self.__lonDimension) and variable in lonVariables:
                variables.append(variable)

        self.__variables = variables

        return variables

    def __is0to360(self):
        '''
        We need to check to get the Min and Max lon values to determine if
        the dataset goes from 0 - 360 or -180 - 180
        '''

        is0to360 = False

        netCDFProp = self.__ncFileProperties;
        dimLen = netCDFProp.getDimensionSize(self.__lonDimension);

        # Determining the min and max dim.  Looking for values less than 0 and/or greater
        # than 180
        minLon = 360
        maxLon = -180
        start = time.time()

        for index in range(0, dimLen):

            value = netCDFProp.getDimensionValue(self.__lonDimension, index)
            elapsed = (time.time() - start)
            arcpy.AddMessage("Get Dimension Value for index: " + str(index) + "...." + str(elapsed))

            if value < minLon:
                minLon = value;
            if value > maxLon:
                maxLon = value;

        arcpy.AddMessage("Min Lon Value: " + str(minLon))
        arcpy.AddMessage("Max Lon Value: " + str(maxLon))

        # Checking if the min and or max values fall in the right range
        if minLon >= 0 and maxLon > 180:
            is0to360 = True

        arcpy.AddMessage("Is 0 to 360" + str(is0to360))

        return is0to360

    def getDimensions(self):
        return self.__dimensions

    def getVariables(self):
        return self.__variables

    def getLatDimension(self):
        return self.__latDimension

    def getLonDimension(self):
        return self.__lonDimension

    def getLatValue(self):
        return self.__latValue

    def getLonValue(self):
        return self.__lonValue

    def makeNetCDFTable(self, inpnt, variableName, rowDim, outTableView):
        '''
         Method Name:  makeNetCDFTable
         Description:  Creates a netCDF Table view from the first point selected
                       on the map, the variable choosen and the row dimension.  The
                       table contains every variable value at that particular point
                       for each slice of the row dimension.  For example: Every sea
                       temperature (Variable) value at each elevation(Row Dimension)
                       for a particular point(Input Point) in sea.
         Input:
                       inpnt:         The selected Point
                       variableName:  The variable to get each value
                       rowDim:        How to slice up the netCDF file
                       outTableView:  The table being outputed
        '''

        netCDFSource = str(self.__sourceLocation)

        lonVar = self.getLonDimension()
        latVar = self.getLatDimension()

        # Printing out the inputs
        arcpy.AddMessage("Input NetCDF: " + netCDFSource)
        arcpy.AddMessage("Variable Name: " + str(variableName))
        arcpy.AddMessage("Row Dim: " + str(rowDim))
        arcpy.AddMessage("Lat Dim: " + latVar)
        arcpy.AddMessage("Lon Dim: " + lonVar)

        copyStartTime = time.time()

        arcpy.CopyFeatures_management(inpnt, r'in_memory\updateFeat')

        elapsedTime = (time.time() - copyStartTime)
        arcpy.AddMessage("Copy Features " + str(elapsedTime))

        is0to360 = self.__is0to360()

        # Only getting the first point.  Others ignored
        with arcpy.da.SearchCursor('in_memory\updateFeat', ('SHAPE@X', 'SHAPE@Y')) as cursor:
            for row in cursor:

                # Store x,y coordinates of current point
                if is0to360:
                    self.__lonValue = row[0] + 180
                else:
                    self.__lonValue = row[0]

                self.__latValue = row[1]

                arcpy.AddMessage(lonVar + ": " + str(self.__lonValue) + " " + latVar + ": " + str(self.__latValue))

        netCDFStartTime = time.time()
        arcpy.MakeNetCDFTableView_md(netCDFSource, variableName, outTableView, rowDim,
                                     lonVar + " " + str(self.__lonValue) + ";" + latVar + " " + str(self.__latValue),
                                     "BY_VALUE")

        elapsedNetCDFTime = (time.time() - copyStartTime)
        arcpy.AddMessage("Make netCDF Table Time: " + str(elapsedNetCDFTime))

    @staticmethod
    def isNetCDF(netCDFLoc):
        '''
         Method Name:  isNetCDF
         Description:  Checks to see if the file is a netCDF by checking that it ends
                       with '.nc'
         Input:        NetCDF File location
         Output:       True/False if the file is a NetCDF File
        '''
        isNetCDFBool = True

        if not str(netCDFLoc).endswith(".nc"):
            isNetCDFBool = False

        return isNetCDFBool

    @staticmethod
    def getNetCDFPathfromLayer(netCDFLayer):
        '''
         Method Name:  getNetCDFPathfromLayer
         Description:  The data source of a NetCDF Raster layer contains some sort of
                       additional text that represents the raster in memory.
                       We want to strip that piece off so that we only have the
                       exact loaction of the netCDF file
         Input:        NetCDF Raster Layer
         Output:       String: The path to the NetCDF File
        '''
        datasource = str(netCDFLayer.dataSource)
        return datasource.replace('\\' + netCDFLayer.datasetName, '')


def raster_to_xyz(raster, raster_name, output, no_data_value):
    raster_desc = arcpy.Describe(raster)
    raster_x_min = raster_desc.Extent.XMin
    raster_y_min = raster_desc.Extent.YMin
    raster_ch = raster_desc.MeanCellHeight
    raster_cw = raster_desc.MeanCellWidth
    raster_h = raster_desc.Height
    raster_w = raster_desc.Width

    # Adjust start of xmin and ymin to adjust for centroids
    raster_x_min_adj = raster_x_min + (0.5 * raster_cw)
    raster_y_min_adj = raster_y_min + (0.5 * raster_ch)

    ## Build coordinates
    def raster_centroid_x(raster_h, raster_w):
        coordinates_x = raster_x_min_adj + (raster_cw * raster_w)
        return coordinates_x
    raster_coordinates_x = np.fromfunction(raster_centroid_x, (raster_h, raster_w))  # numpy array of X coord

    def raster_centroid_y(raster_h, raster_w):
        coordinates_y = raster_y_min_adj + (raster_ch * raster_h)
        return coordinates_y
    raster_coordinates_y = np.fromfunction(raster_centroid_y, (raster_h, raster_w))  # numpy array of Y coord

    coordinates_y = raster_y_min_adj + (raster_ch * raster_h)
    coordinates_x = raster_x_min_adj + (raster_cw * raster_w)

    # combine arrays of coordinates (although array for Y is before X, dstack produces [X, Y] pairs)
    raster_coordinates = np.dstack((raster_coordinates_y, raster_coordinates_x))
    del raster_coordinates_y, raster_coordinates_x
    gc.collect()

    ##Raster conversion to NumPy Array
    raster_values = arcpy.RasterToNumPyArray(raster)

    # flip array upside down - then lower left corner cells has the same index as cells in coordinates array
    raster_values = np.flipud(raster_values)

    # combine coordinates and value
    raster_values_full = np.dstack((raster_coordinates, raster_values))
    out = csv.writer(open(os.path.join(output, str(raster_name) + ".yxz"), "wb"),
                     delimiter=' ', quoting=csv.QUOTE_NONNUMERIC)
    (height, width, dim) = raster_values_full.shape
    for row in range(0, height):
        for col in range(0, width):
            out.writerow(raster_values_full[row, col])

    del raster_values, raster_values_full, raster_coordinates,
    gc.collect()

    return coordinates_x, coordinates_y



def main():
    tool = DeepSeaSDMToolsExtractWOANetCDF_mp()
    tool.execute(tool.getParameterInfo(), None)


if __name__ == '__main__':
    main()


