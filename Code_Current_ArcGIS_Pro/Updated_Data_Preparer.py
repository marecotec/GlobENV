 #!/usr/bin/env python
import arcpy
import os
from arcpy import env
from Includes import load_depth_string
from Includes import NetCDFFile
from Includes import raster_to_xyz
import pandas as pd
import gc
import numpy as np
import time

gc.enable()

arcpy.CheckOutExtension("Spatial")


t_start = time.clock()

arcpy.env.overwriteOutput = True

input_woa_netcdf = r"G:\GlobENV_SRTM15\Processing\1_Environment_Layers\glodapv2\Raw_Data\GLODAPv2.2016b_MappedClimatologies\GLODAPv2.2016b.pHtsinsitutp.nc"
variable_name = 'pHtsinsitutp'
lat_name = "lat"
lon_name = "lon"
depth_name = "depth_surface"
depths = "WOA05"
interpolation_procedure = "Natural Neighbor and IDW"
interpolation_resolution = "0.2"
extraction_extent = "-190 -100 190 100"
temporary_directory = r"G:\GlobENV_SRTM15\Processing\1_Environment_Layers\glodapv2\pHtsinsitutp_temp"
output_directory = r"G:\GlobENV_SRTM15\Processing\1_Environment_Layers\glodapv2\pHtsinsitutp"
coordinate_system = "PROJCS['World_Mercator',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Mercator'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',0.0],PARAMETER['Standard_Parallel_1',0.0],UNIT['Meter',1.0]]"
createxyz = "Only Projected"

if not os.path.exists(output_directory):
    os.makedirs(output_directory)

if not os.path.exists(os.path.join(output_directory, "Projected")):
    os.makedirs(os.path.join(output_directory, "Projected"))

if not os.path.exists(os.path.join(output_directory, "Geographic")):
    os.makedirs(os.path.join(output_directory, "Geographic"))

if not os.path.exists(temporary_directory):
    os.makedirs(temporary_directory)



arcpy.AddMessage("Extracting " + str(input_woa_netcdf) + ".")

# Set environment variables and build other variables for processing
arcpy.env.mask = ""
arcpy.env.workspace = temporary_directory
depth_range = load_depth_string(depths)

# Process goes: 1) Convert depth layer to a point file. 2) Interpolate to selected resolution using your selected
# interpolation procedure, 3) Save that layer back into a layer with the name of the variable, plus the
# actual depth value associated with it, you will end up with a specified direction of n rasters (n = number of
# depth layers.

# First lets give an indication of the magnitude of this analysis
arcpy.AddMessage("There are " + str(len(depth_range)) + " depths to process.")

count_geo = 0
count_proj = 0

depth_count = 0

for i in depth_range:
    arcpy.AddMessage("Working on " + str(int(i)))

    # Set some values that we will use to extract data from the NetCDF file

    dimensionValues = str(depth_name) + " " + str(int(depth_count))

    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(4326)
    arcpy.env.extent = None

    # 1 Extract layer to a temporary feature class
    arcpy.MakeNetCDFFeatureLayer_md(in_netCDF_file=input_woa_netcdf, variable=variable_name, x_variable=str(lon_name),
                                    y_variable=str(lat_name),
                                    out_feature_layer="netcdflayer",
                                    row_dimension=str(lat_name) + ";" + str(lon_name),
                                    z_variable="", m_variable="", dimension_values=dimensionValues,
                                    value_selection_method="BY_VALUE")

    arcpy.CopyFeatures_management("netcdflayer", os.path.join(temporary_directory, "out_feature1_" + str(int(i)) + ".shp"))

    arcpy.Delete_management("netcdflayer")

    out_temp_layer = os.path.join(temporary_directory, "out_feature1_" + str(int(i)) + ".shp")

    depth_count = depth_count + 1

    # as someone decided to use 20-380 as the coordinate ref..
    with arcpy.da.UpdateCursor(out_temp_layer, 'SHAPE@XY') as cursor:
        for row in cursor:
            if row[0][0] > 180:
                 cursor.updateRow([[row[0][0] - 360,
                                   row[0][1]]])

    # 2 Interpolate to higher resolution and 3 save to output director
    arcpy.CopyFeatures_management(out_temp_layer, os.path.join(temporary_directory, "out_feature2_" + str(int(i)) + ".shp"))
    #arcpy.Delete_management(out_temp_layer)
    out_temp_layer = os.path.join(temporary_directory, "out_feature2_" + str(int(i)) + ".shp")
    arcpy.RecalculateFeatureClassExtent_management(out_temp_layer)

    arcpy.env.extent = extraction_extent

    if interpolation_procedure == "IDW":
        arcpy.AddMessage("Interpolating " + str(int(i)) + " using IDW")
        arcpy.gp.Idw_sa(out_temp_layer, variable_name[0:10],
                        os.path.join(output_directory, "Geographic", variable_name[0:4] + str(int(i))),
                        interpolation_resolution, "2", "VARIABLE 10", "")
    elif interpolation_procedure == "Spline":
        arcpy.AddMessage("Interpolating " + str(int(i)) + " using Spline")
        arcpy.gp.Spline_sa(out_temp_layer, variable_name[0:10],
                           os.path.join(output_directory, "Geographic", variable_name[0:4] + str(int(i))),
                           interpolation_resolution, "TENSION", "0.1", "10")
        arcpy.Delete_management(os.path.join(output_directory, "out.shp"))
    elif interpolation_procedure == "Kriging":
        arcpy.AddMessage("Interpolating " + str(int(i)) + " using Ordinary Kriging")
        arcpy.gp.Kriging_sa(out_temp_layer, variable_name[0:10],
                            os.path.join(output_directory, "Geographic", variable_name[0:4] + str(int(i))),
                            "Spherical " + str(interpolation_resolution), interpolation_resolution,
                            "VARIABLE 10", "")

    elif interpolation_procedure == "Natural Neighbor":
        arcpy.AddMessage("Interpolating " + str(int(i)) + " using Natural Neighbor")
        arcpy.NaturalNeighbor_3d(out_temp_layer, variable_name[0:10],
                                 os.path.join(output_directory, "Geographic", variable_name[0:4] + str(int(i))),
                                 interpolation_resolution)

    elif interpolation_procedure == "Natural Neighbor and IDW":
        arcpy.AddMessage("Interpolating " + str(int(i)) + " using Natural Neighbor and IDW")

        if not os.path.exists(os.path.join(output_directory, "temp", "idw")):
            os.makedirs(os.path.join(output_directory, "temp", "idw"))

        if not os.path.exists(os.path.join(output_directory, "temp", "nat")):
            os.makedirs(os.path.join(output_directory, "temp", "nat"))

        arcpy.NaturalNeighbor_3d(out_temp_layer, variable_name[0:10],
                                 os.path.join(output_directory, "temp", "nat", variable_name[0:4] + str(int(i))),
                                 interpolation_resolution)

        arcpy.gp.Idw_sa(out_temp_layer, variable_name[0:10],
                        os.path.join(output_directory, "temp", "idw", variable_name[0:4] + str(int(i))),
                        interpolation_resolution, "2", "VARIABLE 10", "")

        input_rasters = [os.path.join(output_directory, "temp", "nat", variable_name[0:4] + str(int(i))),
                         os.path.join(output_directory, "temp", "idw", variable_name[0:4] + str(int(i)))]

        arcpy.MosaicToNewRaster_management(input_rasters=input_rasters,
                                           output_location=os.path.join(output_directory, "Geographic"),
                                           raster_dataset_name_with_extension=variable_name[0:4] + str(int(i)),
                                           coordinate_system_for_the_raster="",
                                           pixel_type="8_BIT_UNSIGNED", cellsize="",
                                           number_of_bands="1", mosaic_method="FIRST", mosaic_colormap_mode="FIRST")

    elif interpolation_procedure == "None":
        arcpy.AddMessage("Making a raster for " + str(int(i)))
        arcpy.MakeNetCDFRasterLayer_md(in_netCDF_file=input_woa_netcdf, variable=variable_name,
                                       x_dimension=lon_name, y_dimension=lat_name,
                                       out_raster_layer=variable_name[0:4] + str(int(i)),
                                       band_dimension="", dimension_values="",
                                       value_selection_method="BY_VALUE")
        arcpy.CopyRaster_management(variable_name[0:4] + str(int(i)),
                                    os.path.join(output_directory, "Geographic", variable_name[0:4] + str(int(i))),
                                    "", "", "", "NONE", "NONE", "")

    #arcpy.Delete_management(out_temp_layer)

    if len(coordinate_system) > 1:
        arcpy.AddMessage("Reprojecting " + variable_name[0:4] + str(int(i)) + ".")
        arcpy.ProjectRaster_management(os.path.join(output_directory, "Geographic", variable_name[0:4] + str(int(i))),
                                       os.path.join(output_directory, "Projected", variable_name[0:4] + str(int(i))),
                                       coordinate_system, "NEAREST", "#", "#", "#", "#")

    arcpy.AddMessage("Generating master file for trilinear interpolation")

    if createxyz == "Only Geographic" or createxyz == "Both":
        if not os.path.exists(os.path.join(output_directory, "Geographic_yxz")):
            os.makedirs(os.path.join(output_directory, "Geographic_yxz"))
        raster_to_xyz(os.path.join(output_directory, "Geographic", variable_name[0:4] + str(int(i))),
                      variable_name[0:4] + str(int(i)),
                      os.path.join(output_directory, "Geographic_yxz"), 349000000.0)
        depth = int(filter(str.isdigit, str(i)))

        if count_geo == 0:
            df = pd.read_csv(os.path.join(output_directory, "Geographic_yxz", variable_name[0:4] + str(int(i)) + ".yxz"),
                             header=0, names=["y", "x", "z"], sep=" ", dtype={"y": np.float64,  "x": np.float64, "z": np.float64})
            master = df[["x", "y", "z"]].copy()
            master.columns = ["x", "y", int(depth)]
            master.to_pickle(os.path.join(output_directory, "Geographic_yxz", "master.pkl"))
            os.remove(os.path.join(output_directory, "Geographic_yxz", variable_name[0:4] + str(int(i)) + ".yxz"))
            del df, master
            gc.collect()
            count_geo = 1
        if count_geo == 1:
            master = pd.read_pickle(os.path.join(output_directory, "Geographic_yxz", "master.pkl"))
            df = pd.read_csv(os.path.join(output_directory, "Geographic_yxz", variable_name[0:4] + str(int(i)) + ".yxz"),
                             header=0, names=["y", "x", "z"], sep=" ", dtype={"y": np.float64,  "x": np.float64, "z": np.float64})
            master[int(depth)] = df["z"].copy()
            master.to_pickle(os.path.join(output_directory, "Geographic_yxz", "master.pkl"))
            os.remove(os.path.join(output_directory, "Geographic_yxz", variable_name[0:4] + str(int(i)) + ".yxz"))
            del df, master
            gc.collect()

    if createxyz == "Only Projected" or createxyz == "Both":
        if not os.path.exists(os.path.join(output_directory, "Projected_yxz")):
            os.makedirs(os.path.join(output_directory, "Projected_yxz"))
        raster_to_xyz(os.path.join(output_directory, "Projected", variable_name[0:4] + str(int(i))),
                      variable_name[0:4] + str(int(i)),
                      os.path.join(output_directory, "Projected_yxz"), 349000000.0)

        depth = int(filter(str.isdigit, str(i)))

        if count_proj == 0:
            df = pd.read_csv(os.path.join(output_directory, "Projected_yxz", variable_name[0:4] + str(int(i)) + ".yxz"),
                             header=0, names=["y", "x", "z"], sep=" ", dtype={"y": np.float64,  "x": np.float64, "z": np.float64})
            master = df[["x", "y"]].copy()
            master.columns = ["x", "y"]
            master = np.round(master, 4)
            master.to_pickle(os.path.join(output_directory, "Projected_yxz", "xy_coords.pkl"))
            master.to_pickle(os.path.join(output_directory, "Projected", "xy_coords.pkl"))
            del master
            gc.collect()
            master_z = df[["z"]].copy()
            master_z.columns = [int(depth)]
            master_z = np.round(master_z, 4)
            master_z.to_pickle(os.path.join(output_directory, "Projected_yxz", str(int(i)) + ".pkl"))
            os.remove(os.path.join(output_directory, "Projected_yxz", variable_name[0:4] + str(int(i)) + ".yxz"))
            del df, master_z
            gc.collect()
            count_proj = 1
        elif count_proj == 1:
            df = pd.read_csv(os.path.join(output_directory, "Projected_yxz", variable_name[0:4] + str(int(i)) + ".yxz"),
                             header=0, names=["y", "x", "z"], sep=" ", dtype={"y": np.float64,  "x": np.float64, "z": np.float64})
            master_z = df[["z"]].copy()
            master_z.columns = [int(depth)]
            master_z = np.round(master_z, 4)
            master_z.to_pickle(os.path.join(output_directory, "Projected_yxz", str(int(i)) + ".pkl"))
            os.remove(os.path.join(output_directory, "Projected_yxz", variable_name[0:4] + str(int(i)) + ".yxz"))
            del df, master_z
            gc.collect()

arcpy.AddMessage("Making pyramids and statistics for outputs")
arcpy.BuildPyramidsandStatistics_management(in_workspace=os.path.join(output_directory, "Geographic"), include_subdirectories="NONE",
                                            build_pyramids="BUILD_PYRAMIDS",
                                            calculate_statistics="CALCULATE_STATISTICS", BUILD_ON_SOURCE="NONE",
                                            block_field="", estimate_statistics="NONE", x_skip_factor="1",
                                            y_skip_factor="1", ignore_values="", pyramid_level="-1",
                                            SKIP_FIRST="NONE", resample_technique="NEAREST",
                                            compression_type="DEFAULT", compression_quality="75",
                                            skip_existing="SKIP_EXISTING")
if len(coordinate_system) > 1:
    arcpy.BuildPyramidsandStatistics_management(in_workspace=os.path.join(output_directory, "Projected"), include_subdirectories="NONE",
                                                build_pyramids="BUILD_PYRAMIDS",
                                                calculate_statistics="CALCULATE_STATISTICS", BUILD_ON_SOURCE="NONE",
                                                block_field="", estimate_statistics="NONE", x_skip_factor="1",
                                                y_skip_factor="1", ignore_values="", pyramid_level="-1",
                                                SKIP_FIRST="NONE", resample_technique="NEAREST",
                                                compression_type="DEFAULT", compression_quality="75",
                                                skip_existing="SKIP_EXISTING")

arcpy.AddMessage("Script complete in %s seconds." % (time.clock() - t_start))