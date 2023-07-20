##################################################################################################
#   ___ _     _    ___ _  ___   __
#  / __| |___| |__| __| \| \ \ / /
# | (_ | / _ \ '_ \ _|| .` |\ V /
#  \___|_\___/_.__/___|_|\_| \_/
#
# GlobENV - Benthic data layer interpolation
#
# Dr Andrew J Davies - University of Rhode Island, USA - github: marecotec
#
##################################################################################################
# Copyright 2019 Dr Andrew J Davies
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

import sys, os, arcpy, multiprocessing, csv, gc
from scipy.interpolate import RegularGridInterpolator
import numpy as np
from arcpy import env
from arcpy.sa import *
import pandas as pd
from functools import partial
import subprocess


def globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, chunk_mode, test_mode, verbose_mode):

    environment_name = environment_name

    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("Spatial")

    print("GlobENV - Trilinear interpolation of environmental data onto bathymetry")
    print("... Processing " + str(input_bathymetry_folder))

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    env.workspace = input_environment

    # Build list of environmental data
    input_environment_list = arcpy.ListRasters("*", "ALL")
    input_environment_depth = []

    # Split name to extract depth
    for i in input_environment_list:
        def remove_prefix(text, prefix):
            if text.startswith(prefix):
                return text[len(prefix):]
            return text  # or whatever

        input_environment_depth.append(int(remove_prefix(i,environment_name)))

    input_environment_name = ''.join([i for i in str(input_environment_list[0]) if not i.isdigit()])

    input_environment_depth.sort()

    print("... There are " + str(len(input_environment_depth)) + " environmental layers")

    if len(input_environment_depth) < 2:
        print("Error: Too few environmental layers for Trilinear interpolation")
        sys.exit(0)
    else:
        pass

    arcpy.env.outputCoordinateSystem = os.path.join(input_bathymetry_folder, "bath_config")
    arcpy.env.cellSize = os.path.join(input_bathymetry_folder, "bath_config")

    # Check for Projected coordinate system - IMPORTANT - FAIL if not with warning.
    input_bathymetry_desc = arcpy.Describe(os.path.join(input_bathymetry_folder, "bath_config"))
    input_bathymetry_sr = input_bathymetry_desc.spatialReference
    input_bathymetry_cs = input_bathymetry_desc.meanCellHeight

    input_environment0_desc = arcpy.Describe(input_environment_list[0])
    input_environment0_sr = input_environment0_desc.spatialReference
    input_environment0_cs = input_environment0_desc.meanCellHeight

    if input_bathymetry_sr.type == "Projected" and input_environment0_sr.type == "Projected":
        print("... Both depth and environmental data are in valid projection")
        if input_bathymetry_sr.name == input_environment0_sr.name:
            pass
        else:
            print("Error: Depth and environmental data are in different projections, must be identical")
            sys.exit(0)
        pass
    else:
        print("Error: Both depth and environmental data need to be in projected coordinate system")
        sys.exit(0)

    if not os.path.exists(os.path.join(output_directory, "Outputs")):
        os.makedirs(os.path.join(output_directory, "Outputs"))

    env.workspace = input_bathymetry_folder
    input_bathymetry_list = arcpy.ListRasters("sp*", "GRID")

    print("... There are " + str(len(input_bathymetry_list)) + " depth blocks to process.")

    input_bathymetry_list_done = []

    rgi = build_globenv_array(input_environment, environment_name)

    for i in input_bathymetry_list:
        try:
            r = arcpy.Raster(os.path.join(input_bathymetry_folder, str(i)))
        except:
            r = False
            print("... Issue reading raster " + str(i))

        del r
        if not os.path.exists(os.path.join(output_directory, "Outputs", str(i) + ".asc")):
            input_bathymetry_list_done.append(i)

    if len(input_bathymetry_list_done) > 0:
        if test_mode:
            print("... Test mode: Processing in a single thread")
            for i in input_bathymetry_list_done:
                mpprocess(output_directory, input_environment_name, input_bathymetry_folder, verbose_mode, rgi, i)
        else:
            print("... Parallel mode: There are " + str(len(input_bathymetry_list_done)) + " depth left to process.")

            if int(len(input_bathymetry_list_done)) < int(cpu_cores_used):
                pool = multiprocessing.Pool(int(len(input_bathymetry_list_done)))
                print("... Will use " + str(int(len(input_bathymetry_list_done))) + " cores for processing")
            else:
                pool = multiprocessing.Pool(int(cpu_cores_used))
                print("... Will use " + str(cpu_cores_used) + " cores for processing")

            ans = pool.map(partial(mpprocess, output_directory, input_environment_name, input_bathymetry_folder, verbose_mode, rgi), input_bathymetry_list_done)
            pool.close()
            pool.join()

    else:
        print("... All rasters processed")


    if chunk_mode == 'gdal':
        env.workspace = os.path.join(output_directory, "Outputs")
        output_list = arcpy.ListRasters("*", "ALL")
        gdalchunk(output_directory, output_list)

    elif chunk_mode == 'arcpy':

        print("... Mosaicking using arcpy, very slow")

        env.workspace = os.path.join(output_directory, "Outputs")
        output_list = arcpy.ListRasters("*", "ALL")

        if not os.path.exists(os.path.join(output_directory, "Outputs_C")):
            os.makedirs(os.path.join(output_directory, "Outputs_C"))

        if len(output_list) > 200:
            output_list_chunk = [output_list[x:x + 20] for x in range(0, len(output_list), 20)]
        else:
            output_list_chunk = [output_list[x:x + 5] for x in range(0, len(output_list), 5)]

        print("... Building pool to mosaic " + str(len(output_list)) + " rasters, " + "in " + str(
            len(output_list_chunk)) + " chunks.")

        if int(len(output_list_chunk)) < int(cpu_cores_used):
            pool2 = multiprocessing.Pool(int(len(output_list_chunk)))
            print("... Will use " + str(int(len(output_list_chunk))) + " cores for mosaic processing")
        else:
            pool2 = multiprocessing.Pool(int(cpu_cores_used))
            print("... Will use " + str(cpu_cores_used) + " cores for mosaic processing")

        func_mosaic = partial(mpchunk, output_directory, input_bathymetry_cs)
        pool2.map(func_mosaic, output_list_chunk)
        pool2.close()
        pool2.join()

        chunk_list_a = []
        for chunk in output_list_chunk:
            if arcpy.Exists(os.path.join(output_directory, "Outputs_C", str(chunk[0]) + "_f", "a")):
                chunk_list_a.append(
                    os.path.join(output_directory, "Outputs_C", str(chunk[0]) + "_f", "a"))

        env.workspace = os.path.join(output_directory, "Outputs_C")

        print("... Mosaicking " + str(len(chunk_list_a)) + " chunks together.")

        arcpy.MosaicToNewRaster_management(
            input_rasters=chunk_list_a,
            output_location=os.path.join(output_directory),
            raster_dataset_name_with_extension="f", coordinate_system_for_the_raster="",
            pixel_type="32_BIT_FLOAT",
            cellsize=input_bathymetry_cs, number_of_bands="1", mosaic_method="MEAN", mosaic_colormap_mode="MATCH")

        arcpy.DefineProjection_management(in_dataset=os.path.join(output_directory, "f"),
                                          coor_system=os.path.join(input_bathymetry_folder, "bath_config"))

        if arcpy.Exists(os.path.join(output_directory, "tri_output")):
            arcpy.Delete_management(os.path.join(output_directory, "tri_output"))

        arcpy.Rename_management(os.path.join(output_directory, "f"),
                                os.path.join(output_directory, "tri_output"))

        print("... GlobENV complete!")

def build_globenv_array(input_environment, environment_name):
    # 1 Extract list of files.
    env_files = []
    env_values = []
    env.workspace = input_environment
    rasterlist1 = arcpy.ListRasters("*")
    for f in rasterlist1:
        f2 = float(f.replace(environment_name, ""))
        env_values.append(f2)
        env_files.append(f)
    else:
        pass

    env_file_list = sorted(zip(env_files, env_values), key=lambda tup: tup[1])

    # 2 Get xyz values needed to build array
    xy_coords = pd.read_pickle(os.path.join(input_environment, "xy_coords.pkl"))

    x_vals = np.unique(xy_coords["x"])
    y_vals = np.unique(xy_coords["y"])

    env_name = []
    env_depth = []

    for pair in env_file_list:
        env_name_i, env_depth_i = pair
        env_name.append(env_name_i)
        env_depth.append(env_depth_i)

    def pad(data):
        bad_indexes = np.isnan(data)
        good_indexes = np.logical_not(bad_indexes)
        good_data = data[good_indexes]
        interpolated = np.interp(bad_indexes.nonzero()[0], good_indexes.nonzero()[0], good_data)
        data[bad_indexes] = interpolated
        return data

    data = np.array([np.flipud(arcpy.RasterToNumPyArray(os.path.join(input_environment, environment_name + "%f" % f),
                                                            arcpy.Point(float(min(x_vals)), float(min(y_vals))), len(x_vals),
                                                            len(y_vals), nodata_to_value=np.nan)) for f in env_depth])
    data = data.T

    z_vals = np.unique(np.asarray(env_depth[::-1], dtype=float))

    print("Environment depth array: " + str(data.size))
    print("Environment depth array (z): " + str(z_vals.size))
    print("Environment depth array (x): " + str(x_vals.size))
    print("Environment depth array (y): " + str(y_vals.size))

    rgi = RegularGridInterpolator((x_vals, y_vals, z_vals), data, method='linear')

    print("RGI Env Grid Created")
    return rgi

def mpprocess(output_directory, input_environment_name, input_bathymetry_folder, verbose_mode, rgi, input_bathymetry_list):
    try:
        # build single environmental grid
        arcpy.CheckOutExtension("Spatial")
        arcpy.env.overwriteOutput = True

        if not os.path.exists(os.path.join(output_directory, "1_Temp", str(input_bathymetry_list))):
            os.makedirs(os.path.join(output_directory, "1_Temp", str(input_bathymetry_list)))

        env.workspace = os.path.join(output_directory, "1_Temp", str(input_bathymetry_list))
        env.scratchWorkspace = os.path.join(output_directory, "1_Temp", str(input_bathymetry_list))

        log_file("GlobENV Logfile\n\n", output_directory, input_bathymetry_list, verbose_mode)
        log_file("Bathymetry: " + str(input_bathymetry_list) + "\n", output_directory, input_bathymetry_list, verbose_mode)
        log_file("Environment name: " + str(input_environment_name) + "\n", output_directory, input_bathymetry_list, verbose_mode)

        try:
            arcpy.RasterToASCII_conversion(os.path.join(input_bathymetry_folder, input_bathymetry_list),
                                           os.path.join(output_directory, str(input_bathymetry_list) + ".asc"))
            input_bathymetry_split = arcpy.Raster(os.path.join(output_directory,str(input_bathymetry_list) + ".asc"))
            input_bathymetry_split_i = input_bathymetry_list
        except:
            log_file("Unable to read bathymetric layer" + "\n", output_directory, input_bathymetry_list, verbose_mode)
            input_bathymetry_split = False
            input_bathymetry_split_i = False

        if arcpy.Exists(input_bathymetry_split) and not \
                os.path.exists(
                    os.path.join(output_directory, "Outputs", str(input_bathymetry_split) + ".asc")):

            no_data_value = 349000000.0
            input_bathymetry_split_2 = Con(IsNull(input_bathymetry_split), no_data_value,
                                           input_bathymetry_split)

            raster_to_xyz(input_bathymetry_split_2, str(input_bathymetry_list),
                          os.path.join(output_directory, "1_Temp", str(input_bathymetry_list)), no_data_value)

            depth_array = pd.read_csv(
                os.path.join(output_directory, "1_Temp", str(input_bathymetry_list), str(input_bathymetry_list) + ".yxz"),
                header=None,
                names=["y", "x", "depth"], sep=" ", low_memory=True)

            depth_array.loc[depth_array["depth"] == no_data_value, "depth"] = np.nan

            depth_array_min = depth_array['depth'].min()

            tri_value_list = []

            if depth_array_min < 0:
                depth_array["depth"] *= -1

            depth_array = depth_array.iloc[::-1]

            log_file("Depth array: " + str(depth_array) + "\n", output_directory,
                     input_bathymetry_list, verbose_mode)

            log_file("Trilinearly interpolating: " + str(input_bathymetry_split_i) + "\n", output_directory,
                     input_bathymetry_list, verbose_mode)

            for i in depth_array.index:
                depth_val = depth_array._get_value(i, "depth")
                if not np.isnan(depth_val):

                        tri_value_list.append(rgi((depth_array._get_value(i, "x"), depth_array._get_value(i, "y"), depth_val)))

                elif np.isnan(depth_val):
                    tri_value_list.append("-9999")

            # for index, row in depth_array.iterrows():
            #     if not np.isnan(row["depth"]):
            #         tri_value_list.append(rgi((row["x"], row["y"], row["depth"])))
            #     elif np.isnan(row["depth"]):
            #         tri_value_list.append("-9999")

            # Get headers
            with open(os.path.join(output_directory,
                                   str(input_bathymetry_split_i) + ".asc")) as myfile:
                header_ascii = myfile.readlines()[0:6]  # put here the interval you want

            # Write headers to new file
            with open(os.path.join(output_directory, "Outputs",
                                   str(input_bathymetry_split_i) + ".asc"), 'a') as ascii_file:
                for wr in header_ascii:
                    ascii_file.write(wr)

            with open(os.path.join(output_directory, "Outputs",
                                   str(input_bathymetry_split_i) + ".asc"), 'a') as ascii_file:
                for tri in tri_value_list:
                    ascii_file.write(str(tri) + " ")
            log_file("Completed block: " + str(input_bathymetry_split_i) + "\n", output_directory,
                     input_bathymetry_list, verbose_mode)

            del tri_value_list, tri, depth_array, depth_array_min, input_bathymetry_split
        else:
            del input_bathymetry_split

    except Exception as e:
        log_file("Recorded error if any: " + str(e) + "\n", output_directory,
                 input_bathymetry_list, verbose_mode)
        log_file("Recorded error if any: " + str(e) + "\n", output_directory,
                 input_bathymetry_list + "_Failed", verbose_mode)
        os.kill(0)

# Stage 3 - Mosaic the chunks back into the extent of the original bathymetry

def gdalchunk(output_directory, output_list_chunk):
    print("... Mosaicking using GDAL")
    arcpy.env.overwriteOutput = True
    if not arcpy.Exists(os.path.join(output_directory, "Outputs_GDAL")):
        os.makedirs(os.path.join(output_directory, "Outputs_GDAL"))
    counter = 0
    for raster_output in output_list_chunk:
        counter += 1
        print(os.path.join(output_directory, "Outputs", raster_output) + " output: " + os.path.join(output_directory, "Outputs_GDAL", "sp" + str(counter) + ".tif"))
        arcpy.Mirror_management(in_raster=os.path.join(output_directory, "Outputs", raster_output),
                                out_raster=os.path.join(output_directory, "Outputs_GDAL", "sp" + str(counter) + ".tif"))

    print("... Calling GDAL")
    gdal_location = r"C:\OSGeo4W\bin\gdalwarp.exe --config GDAL_DATA C:\OSGeo4W\share\gdal "
    gdal_files = os.path.join(output_directory, "Outputs_GDAL", "sp*.tif")
    gdal_output = os.path.join(output_directory, "tri_output.tif")
    with open(os.devnull, "w") as f:
        subprocess.call(gdal_location + gdal_files + " " + gdal_output, stdout=f)
    arcpy.AddMessage("... Mosaicking using GDAL complete")

def mpchunk(output_directory, input_bathymetry_cs, output_list_chunk):
    print("... Processing chunk: " + str(output_list_chunk[0]))
    arcpy.env.overwriteOutput = True
    if not arcpy.Exists(
            os.path.join(output_directory, "Outputs_C", str(output_list_chunk[0]) + "_f", "a")):
        if not os.path.exists(
                os.path.join(output_directory, "Outputs_C", str(output_list_chunk[0]) + "_f")):
            os.makedirs(
                os.path.join(output_directory, "Outputs_C", str(output_list_chunk[0]) + "_f"))

        if not os.path.exists(
                os.path.join(output_directory, "Outputs_C", str(output_list_chunk[0]) + "_f",
                             "temp")):
            os.makedirs(os.path.join(output_directory, "Outputs_C", str(output_list_chunk[0]) + "_f",
                                     "temp"))
        env.scratchWorkspace = os.path.join(output_directory, "Outputs_C",
                                            str(output_list_chunk[0]) + "_f")
        env.workspace = os.path.join(output_directory, "Outputs_C", str(output_list_chunk[0]) + "_f")

        counter = 1
        chunk_list = []
        for raster_output in output_list_chunk:
            if not arcpy.Exists(
                    os.path.join(output_directory, "Outputs_C", str(output_list_chunk[0]) + "_f",
                                 "temp", str(counter))):
                arcpy.Mirror_management(
                    in_raster=os.path.join(output_directory, "Outputs", raster_output),
                    out_raster=os.path.join(output_directory, "Outputs_C",
                                            str(output_list_chunk[0]) + "_f", "temp", str(counter)))
                chunk_list.append(
                    os.path.join(output_directory, "Outputs_C", str(output_list_chunk[0]) + "_f",
                                 "temp", str(counter)))
            counter += 1

        arcpy.MosaicToNewRaster_management(input_rasters=chunk_list,
                                           output_location=os.path.join(output_directory,
                                                                        "Outputs_C",
                                                                        str(output_list_chunk[0]) + "_f"),
                                           raster_dataset_name_with_extension="a",
                                           coordinate_system_for_the_raster="",
                                           pixel_type="32_BIT_FLOAT",
                                           cellsize=input_bathymetry_cs, number_of_bands="1", mosaic_method="MEAN",
                                           mosaic_colormap_mode="MATCH")

    print("... Completed mosaicking chunk: " + str(output_list_chunk[0]))

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

    out = csv.writer(open(os.path.join(output, str(raster_name) + ".yxz"), "w"),
                     delimiter=' ', quoting=csv.QUOTE_NONNUMERIC)
    (height, width, dim) = raster_values_full.shape
    for row in range(0, height):
        for col in range(0, width):
            out.writerow(raster_values_full[row, col])

    del raster_values, raster_values_full, raster_coordinates,
    gc.collect()

    return coordinates_x, coordinates_y

def log_file(text, output_directory, filename, verbose_mode):
    log_file_write = open(os.path.join(output_directory, "Outputs", str(filename) + ".log"), "a")
    log_file_write.write(str(text))
    log_file_write.close()
    if verbose_mode:
        print(str(text))

if __name__ == '__main__':

    # WOA Seasonal
    env_data = [[r"G:\GlobENV_SRTM15\Processing\1_Environment_Layers\steinacher\oc_new_layers", "ocm",
                 "oc"]]

    bathymetries = [
        [r"G:\GlobENV_SRTM15\Processing\2_Bathymetries_and_Outputs\srtm15\Bathymetry\srtm15p_sr", "srtm15"]]

    output_directory_base = r"G:\GlobENV_SRTM15\Processing\2_Bathymetries_and_Outputs\srtm15\stein3oc"

    cpu_cores_used = "28"
    chunk_mode = 'gdal'  # gdal or arcpy - GDAL way faster, need to install GDAL binaries in standard location
    test_mode = False
    verbose_mode = False

    for b in bathymetries:
        for e in env_data:
            if not os.path.exists(os.path.join(output_directory_base, b[1], e[2] + ".tif")):
                globenv(b[0], e[0], e[1], os.path.join(output_directory_base, b[1], e[2]), cpu_cores_used, chunk_mode,
                        test_mode,
                        verbose_mode)
                if chunk_mode == 'gdal':
                    arcpy.Copy_management(os.path.join(output_directory_base, b[1], e[2], "tri_output.tif"),
                                          os.path.join(output_directory_base, b[1], e[2] + ".tif"))
            else:
                print("Skipping: " + b[1] + " " + e[2])