#  GlobENV Data Check
# This compares your bathymetry with your env layers to ensure that GlobENV will run. It checks bottom limits and creates
# dummy data to allow the process to run, it will provide caveats

import arcpy
import os, csv
from arcpy import env

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

def load_environmental_layers(input_env_directory):

    log_file.write("Environmental Data Assessment\n")

    env.workspace = input_env_directory
    # Build list of environmental data
    input_environment_list = arcpy.ListRasters("*", "ALL")
    for i_s in input_environment_list:
        arcpy.management.BuildPyramidsandStatistics(i_s)
    input_environment_depth = []
    # Split name to extract depth
    for i in input_environment_list:
        input_environment_depth.append(int(filter(str.isdigit, str(i))))
    input_environment_name = ''.join([i for i in str(input_environment_list[0]) if not i.isdigit()])
    input_environment_depth.sort()
    log_file.write("Number of environmental data layers: " + str(len(input_environment_depth)) + "\n")
    log_file.write("Environment data name: " + str(input_environment_name) + "\n")
    log_file.write("Environmental depth layers: " + ",".join(str(i) for i in input_environment_depth))
    if len(input_environment_depth) < 2:
        print("Error: Too few environmental layers for Trilinear interpolation. Exiting")
        log_file.write("Error: Too few environmental layers for Trilinear interpolation")
        sys.exit(0)
    else:
        log_file.write("\n\nDescribing first environmental layer properties\n")
        arcpy.CalculateStatistics_management(os.path.join(input_env_directory, input_environment_list[0]))
        env_layer_prop = data_desc_raster(os.path.join(input_env_directory, input_environment_list[0]))

    return env_layer_prop, input_environment_depth, input_environment_name

def load_bathymetry_layer(input_bathymetry):

    log_file.write("\n\nDescribing bathymetry layer properties\n")
    bathy_layer_prop = data_desc_raster(input_bathymetry)
    return bathy_layer_prop

def data_desc_raster(input_layer):

    if input_layer is not None:
        input_layer_prop = {}
        log_file.write("Raster file: " + str(input_layer) + "\n")
        print(str(input_layer))
        arcpy.management.BuildPyramidsandStatistics(input_layer)
        input_layer_desc = arcpy.Describe(input_layer)
        if input_layer_desc.dataType == 'RasterDataset':
            input_layer_prop["extent"] = input_layer_desc.extent
            input_layer_prop["format"] = input_layer_desc.format
            input_layer_prop["cs type"] = input_layer_desc.spatialReference.type
            input_layer_prop["cs name"] = input_layer_desc.spatialReference.name
            input_layer_prop["xy units"] = input_layer_desc.spatialReference.linearUnitName
            input_layer_prop["cell size x"] = arcpy.GetRasterProperties_management(input_layer, "CELLSIZEX").getOutput(0)
            input_layer_prop["cell size y"] = arcpy.GetRasterProperties_management(input_layer, "CELLSIZEY").getOutput(0)


            input_layer_prop["minimum"] = arcpy.GetRasterProperties_management(input_layer, "MINIMUM").getOutput(0)
            input_layer_prop["maximum"] = arcpy.GetRasterProperties_management(input_layer, "MAXIMUM").getOutput(0)
            input_layer_prop["nrows"] = arcpy.GetRasterProperties_management(input_layer, "ROWCOUNT").getOutput(0)
            input_layer_prop["ncols"] = arcpy.GetRasterProperties_management(input_layer, "COLUMNCOUNT").getOutput(0)

            if input_layer_prop.get("cell size x", "1") != input_layer_prop.get("cell size y", "2"):
                print("Warning: Bathymetry x and y cell sizes should be identical (i.e. data must be in square cells).")
                log_file.write("Warning: Bathymetry x and y cell sizes should be identical (i.e. data must be in square cells).\n")

            log_file.write("Raster format: " + str(input_layer_desc.format) + "\n")
            log_file.write("Raster CS type: " + str(input_layer_desc.spatialReference.type) + "\n")
            log_file.write("Raster CS name: " + str(input_layer_desc.spatialReference.name) + "\n")
            log_file.write("Raster xy units: " + str(input_layer_desc.spatialReference.linearUnitName) + "\n")
            log_file.write("Raster cell size: " + str(round(float(input_layer_prop.get("cell size x", "!!!")), 2)) + "\n")
            log_file.write("Raster rows: " + str(input_layer_prop.get("nrows", "!!!")) + "\n")
            log_file.write("Raster cols: " + str(input_layer_prop.get("ncols", "!!!")) + "\n")
            log_file.write("Raster estimated cell count: " + str(int(input_layer_prop.get("nrows", "!!!")) * int(input_layer_prop.get("ncols", "!!!"))) + "\n")
            log_file.write("Value minimum: " + str(float(input_layer_prop.get("minimum", "!!!"))) + "\n")
            log_file.write("Value maximum: " + str(float(input_layer_prop.get("maximum", "!!!"))) + "\n")
        else:
            log_file.write(r"Error: data is not in raster format, not describing.\n")
            print(r"Error: data is not in raster format, not describing.\n")
    else:
        print("Error: file read unsuccessful.")

    return input_layer_prop

def globenv_coordinatesys_check(env_layer_prop, bathy_layer_prop, input_env_directory, temp_dir):

    log_file.write("\nCheck - Coordinate System\n")
    if env_layer_prop.get("cs name", "!!!") == bathy_layer_prop.get("cs name", "!!!"):
        log_file.write("Coordinate system match: True\n")
    else:
        log_file.write("Coordinate system match: False\n")
        print("Error: Coordinate systems of environmental and bathymetry data do not match.")

    log_file.write("\nCheck - Extent\n")
    extent_overlap = env_layer_prop.get("extent", "!!!").overlaps(bathy_layer_prop.get("extent", "!!!"))
    if extent_overlap == True:
        log_file.write("Extents overlap: True\n")
    else:
        log_file.write("Extents overlap: False\n")
        print("Error: Extents do not overlap. Taking remedial action")

        env.workspace = input_env_directory
        input_environment_list = arcpy.ListRasters("*", "ALL")

        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)

        for i in input_environment_list:
            cell_size = float(env_layer_prop.get("cell size x", "!!!"))

            arcpy.management.Shift(i, os.path.join(temp_dir, i + "_ru"), cell_size * 3, cell_size*3, i)
            arcpy.management.Shift(i, os.path.join(temp_dir, i + "_rd"), cell_size * 3, cell_size * -3, i)
            arcpy.management.Shift(i, os.path.join(temp_dir, i + "_ld"), cell_size * -3, cell_size * -3, i)
            arcpy.management.Shift(i, os.path.join(temp_dir, i + "_lu"), cell_size * -3, cell_size * 3, i)

            arcpy.MosaicToNewRaster_management([os.path.join(temp_dir, i + "_ru"), os.path.join(temp_dir, i + "_rd"),
                                                os.path.join(temp_dir, i + "_ld"), os.path.join(temp_dir, i + "_lu"),
                                                i], temp_dir, i, "", "", "", 1, "LAST", "LAST")

            arcpy.Delete_management(os.path.join(temp_dir, i + "_ru"))
            arcpy.Delete_management(os.path.join(temp_dir, i + "_rd"))
            arcpy.Delete_management(os.path.join(temp_dir, i + "_ld"))
            arcpy.Delete_management(os.path.join(temp_dir, i + "_lu"))



    return

def globenv_depthrange_check(env_layer_prop, bathy_layer_prop, input_environment_depth, input_environment_name, input_env_directory):

    log_file.write("\nCheck - Depth ranges\n")

    depth_min = float(bathy_layer_prop.get("minimum", "!!!"))
    depth_max = float(bathy_layer_prop.get("maximum", "!!!"))
    if depth_min < 0.:
        depth_min *= -1
    if depth_max < 0.:
        depth_max *= -1

    if depth_min > depth_max:
        depth_min, depth_max = depth_max, depth_min

    env_depth_min = min(input_environment_depth)
    env_depth_max = max(input_environment_depth)

    depth_min = int(round(depth_min, 0))
    depth_max = int(round(depth_max, 0))

    #env depth should be deeper and shallower than the bathy depth
    if env_depth_min < depth_min:
        log_file.write("Minimum depths: Environment depth is lower than bathymetry (" + str(env_depth_min) + "/" + str(depth_min) + "), no action needed\n")
    elif env_depth_min == depth_min:
        log_file.write("Minimum depths: Environment depth is equal to bathymetry (" + str(env_depth_min) + "/" + str(depth_min) + "), no action needed\n")
    elif env_depth_min > depth_min:
        log_file.write("Minimum depths: Environment depth is greater than bathymetry (" + str(env_depth_min) + "/" + str(depth_min) + "), remedial action taken\n")
        if depth_min == 0:
            arcpy.CopyRaster_management(os.path.join(input_env_directory, input_environment_name + str(env_depth_min)), os.path.join(input_env_directory, input_environment_name + "0"))
            log_file.write("Written dummy env layer: " + str(input_environment_name + "0\n"))
        else:

            arcpy.CopyRaster_management(os.path.join(input_env_directory, input_environment_name + str(env_depth_min)), os.path.join(input_env_directory, input_environment_name + str(depth_min - 1)))
            log_file.write("Written dummy env layer: " + str(round(depth_min - 1)) + "\n")

    if env_depth_max > depth_max:
        log_file.write("Maximum depths: Environment depth is deeper than bathymetry (" + str(env_depth_max) + "/" + str(
            depth_max) + "), no action needed\n")
    elif env_depth_max == depth_max:
        log_file.write("Maximum depths: Environment depth is equal to bathymetry (" + str(env_depth_max) + "/" + str(
            depth_max) + "), no action needed\n")
    elif env_depth_max < depth_max:
        log_file.write(
            "Maximum depths: Environment depth is shallower than bathymetry (" + str(env_depth_min) + "/" + str(
                depth_min) + "), remedial action taken\n")
        depth_difference_min = depth_min - env_depth_min
        depth_difference_max = depth_max - env_depth_max
        log_file.write("Depth difference (min): " + str(depth_difference_min) + "\n")
        log_file.write("Depth difference (max): " + str(depth_difference_max) + "\n")

        arcpy.CopyRaster_management(os.path.join(input_env_directory, input_environment_name + str(env_depth_max)),
                              os.path.join(input_env_directory, input_environment_name + str(env_depth_max + 1)))
        arcpy.CopyRaster_management(os.path.join(input_env_directory, input_environment_name + str(env_depth_max)),
                              os.path.join(input_env_directory, input_environment_name + str(int(round(depth_max)) + 100)))
        log_file.write("Written dummy env layer: " + str(input_environment_name + str(env_depth_max + 1)) + "\n")
        log_file.write("Written dummy env layer: " + str(input_environment_name + str(int(round(depth_max)) + 100)) + "\n")

    return

def generate_xycoords(input_env_directory, ):

    env.workspace = input_env_directory
    input_environment_list = arcpy.ListRasters("*", "ALL")

    def raster_to_xyz(raster, input_env_directory):
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
        out = csv.writer(open(os.path.join(input_env_directory, "xycoords.yxz"), "w"),
                         delimiter=' ', quoting=csv.QUOTE_NONNUMERIC)
        (height, width, dim) = raster_values_full.shape

        for row in range(0, height):
            for col in range(0, width):
                out.writerow(raster_values_full[row, col])

        del raster_values, raster_values_full, raster_coordinates,
        gc.collect()

        return coordinates_x, coordinates_y

    raster_to_xyz(input_environment_list[0], input_env_directory)

    df = pd.read_csv(os.path.join(os.path.join(input_env_directory, "xycoords.yxz")),
                     header=0, names=["y", "x", "z"], sep=" ",
                     dtype={"y": np.float32, "x": np.float32, "z": np.float32})
    master = df[["x", "y"]].copy()
    master.columns = ["x", "y"]
    master = np.round(master, 4)
    master.to_pickle(os.path.join(input_env_directory, "xy_coords.pkl"))
    del master, df
    gc.collect()


def globenv_check(input_env_directory, input_bathymetry, temp_dir):
    print("GlobENV: Running checks on your input environmental data and bathymetry..")
    print(".. check your resulting log file for details, only errors printed here.")
    env_layer_prop, input_environment_depth, input_environment_name = load_environmental_layers(input_env_directory)
    bathy_layer_prop = load_bathymetry_layer(input_bathymetry)
    globenv_depthrange_check(env_layer_prop, bathy_layer_prop, input_environment_depth, input_environment_name, input_env_directory)
    globenv_coordinatesys_check(env_layer_prop, bathy_layer_prop, input_env_directory, temp_dir)
    generate_xycoords(temp_dir)
    return

if __name__ == '__main__':

    verbose_mode = True

    #Input List = input_env_directory, input_bathymetry, report_file

    input_list = [
        #SRTM
        [r"H:\GlobENV_SRTM15\Processing\1_Environment_Layers\world-ocean-atlas-2018\aoxu\Output\Projected",
         r"H:\GlobENV_SRTM15\Processing\2_Bathymetries_and_Outputs\gebco2022\gebco_2022_ocean_behr.tif",
         r"H:\GlobENV_SRTM15\Processing\1_Environment_Layers\world-ocean-atlas-2018\aoxu\Output\report.txt",
         r"H:\GlobENV_SRTM15\Processing\1_Environment_Layers\world-ocean-atlas-2018\aoxu\Output\Projected_Behr"]
        ]


    for i in input_list:
        log_file = open(i[2], "a")
        log_file.write("GlobENV Data Check Logfile\n\n")
        log_file.write("Environmental Data Folder = " + str(i[0]) + "\nBathymetry Dataset = " + str(i[1]) + "\n\n")
        globenv_check(i[0], i[1], i[3])
