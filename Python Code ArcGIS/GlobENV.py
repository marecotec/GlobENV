##################################################################################################
# GlobENV - GlobENV.py
#
# Dr Andy Davies - Bangor University, Wales - github: marecotec
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

import sys, os, arcpy, time, multiprocessing, csv, gc
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import time
from arcpy import env
from arcpy.sa import *
import pandas as pd
from functools import partial


def globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode, verbose_mode):

    verbose_mode = verbose_mode
    environment_name = environment_name

    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("Spatial")

    arcpy.AddMessage("GlobENV - Trilinear interpolation of environmental data onto bathymetry")
    arcpy.AddMessage("... processing " + str(input_bathymetry_folder))

    t_start = time.clock()

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    env.workspace = input_environment

    # Build list of environmental data
    input_environment_list = arcpy.ListRasters("*", "ALL")
    input_environment_depth = []

    # Split name to extract depth
    for i in input_environment_list:
        input_environment_depth.append(int(filter(str.isdigit, str(i))))

    input_environment_name = ''.join([i for i in str(input_environment_list[0]) if not i.isdigit()])

    input_environment_depth.sort()

    arcpy.AddMessage("There are " + str(len(input_environment_depth)) + " environmental layers")

    if len(input_environment_depth) < 2:
        arcpy.AddMessage("Error: Too few environmental layers for Trilinear interpolation")
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
    input_environment0_cs_x_min = input_environment0_desc.extent.XMin + (0.5 * input_environment0_cs)
    input_environment0_cs_y_min = input_environment0_desc.extent.YMin + (0.5 * input_environment0_cs)
    input_environment0_cs_x_max = input_environment0_desc.extent.XMax - (0.5 * input_environment0_cs)
    input_environment0_cs_y_max = input_environment0_desc.extent.YMax - (0.5 * input_environment0_cs)

    if input_bathymetry_sr.type == "Projected" and input_environment0_sr.type == "Projected":
        arcpy.AddMessage("Both depth and environmental data are in valid projection")
        if input_bathymetry_sr.name == input_environment0_sr.name:
            pass
        else:
            arcpy.AddMessage(
                "Error: Depth and environmental data are in different projections, must be identical")
            sys.exit(0)
        pass
    else:
        arcpy.AddMessage("Error: Both depth and environmental data need to be in projected coordinate system")
        sys.exit(0)
    # del input_bathymetry_desc, input_bathymetry_sr, input_environment0_desc, input_environment0_sr

    if not os.path.exists(os.path.join(output_directory, "Outputs")):
        os.makedirs(os.path.join(output_directory, "Outputs"))

    env.workspace = input_bathymetry_folder
    input_bathymetry_list = arcpy.ListRasters("sp*", "GRID")

    arcpy.AddMessage("There are " + str(len(input_bathymetry_list)) + " depth blocks to process.")

    input_bathymetry_list_done = []

    for i in input_bathymetry_list:
        try:
            r = arcpy.Raster(os.path.join(input_bathymetry_folder, str(i)))
        except:
            r = False
            arcpy.AddMessage("Issue reading raster " + str(i))

        del r
        if not os.path.exists(os.path.join(output_directory, "Outputs", str(i) + ".asc")):
            input_bathymetry_list_done.append(i)

    if len(input_bathymetry_list_done) > 0:
        if test_mode:
            arcpy.AddMessage("Test mode: Processing in a single thread")
            for i in input_bathymetry_list_done:
                mpprocess(output_directory, input_environment_depth, input_environment0_cs, input_environment,
                           input_environment_name, input_environment0_cs_x_min, input_environment0_cs_y_min,
                           input_environment0_cs_x_max, input_environment0_cs_y_max, input_bathymetry_folder, test_mode, verbose_mode, i)
        else:
            arcpy.AddMessage("Parallel mode: There are " + str(len(input_bathymetry_list_done)) + " depth left to process.")

            if int(len(input_bathymetry_list_done)) < int(cpu_cores_used):
                pool = multiprocessing.Pool(int(len(input_bathymetry_list_done)))
                arcpy.AddMessage("Will use " + str(int(len(input_bathymetry_list_done))) + " cores for processing")
            else:
                pool = multiprocessing.Pool(int(cpu_cores_used))
                arcpy.AddMessage("Will use " + str(cpu_cores_used) + " cores for processing")

            func = 1
            ans = pool.map(partial(mpprocess, output_directory, input_environment_depth, input_environment0_cs, input_environment,
                           input_environment_name, input_environment0_cs_x_min, input_environment0_cs_y_min,
                           input_environment0_cs_x_max,
                           input_environment0_cs_y_max, input_bathymetry_folder, test_mode, verbose_mode), input_bathymetry_list_done)
            pool.close()
            pool.join()

    else:
        arcpy.AddMessage("All rasters processed")


    if test_mode:
        arcpy.AddMessage("Test mode complete in %s minutes." % ((time.clock() - t_start) / 60.0))
        sys.exit(0)
    else:
        env.workspace = os.path.join(output_directory, "Outputs")
        output_list = arcpy.ListRasters("*", "ALL")

        if not os.path.exists(os.path.join(output_directory, "Outputs_C")):
            os.makedirs(os.path.join(output_directory, "Outputs_C"))

        if len(output_list) > 200:
            output_list_chunk = [output_list[x:x + 20] for x in xrange(0, len(output_list), 20)]
        else:
            output_list_chunk = [output_list[x:x + 5] for x in xrange(0, len(output_list), 5)]

        arcpy.AddMessage("Building pool to mosaic " + str(len(output_list)) + " rasters, " + "in " + str(
            len(output_list_chunk)) + " chunks.")

        if int(len(output_list_chunk)) < int(cpu_cores_used):
            pool2 = multiprocessing.Pool(int(len(output_list_chunk)))
            arcpy.AddMessage("Will use " + str(int(len(output_list_chunk))) + " cores for processing")
        else:
            pool2 = multiprocessing.Pool(int(cpu_cores_used))
            arcpy.AddMessage("Will use " + str(cpu_cores_used) + " cores for processing")

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

        arcpy.AddMessage("Mosaicking " + str(len(chunk_list_a)) + " chunks together.")

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

        arcpy.AddMessage("Output file is called *tri_output*..")
        arcpy.AddMessage("Script complete in %s minutes." % ((time.clock() - t_start) / 60.0))


def mpprocess(output_directory, input_environment_depth, input_environment0_cs, input_environment,
                       input_environment_name, input_environment0_cs_x_min, input_environment0_cs_y_min,
                       input_environment0_cs_x_max,
                       input_environment0_cs_y_max, input_bathymetry_folder, test_mode, verbose_mode, input_bathymetry_list):

    arcpy.CheckOutExtension("Spatial")
    arcpy.env.overwriteOutput = True
    depth_array_min_comp = -9999
    depth_array_max_comp = -9999

    if not os.path.exists(os.path.join(output_directory, "1_Temp", str(input_bathymetry_list))):
        os.makedirs(os.path.join(output_directory, "1_Temp", str(input_bathymetry_list)))

    env.workspace = os.path.join(output_directory, "1_Temp", str(input_bathymetry_list))
    env.scratchWorkspace = os.path.join(output_directory, "1_Temp", str(input_bathymetry_list))

    try:
        arcpy.AddMessage("Reading bathymetric layer: " + str(input_bathymetry_list))
        arcpy.RasterToASCII_conversion(os.path.join(input_bathymetry_folder, input_bathymetry_list),
                                       os.path.join(output_directory,
                                                    str(
                                                        input_bathymetry_list) + ".asc"))
        input_bathymetry_split = arcpy.Raster(os.path.join(output_directory,
                                                           str(
                                                               input_bathymetry_list) + ".asc"))
        input_bathymetry_split_i = input_bathymetry_list
    except:
        arcpy.AddWarning("Unable to read bathymetric layer: " + str(input_bathymetry_list))
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

        input_bathymetry_extent = input_bathymetry_split_2.extent

        depth_array = pd.read_csv(
            os.path.join(output_directory, "1_Temp", str(input_bathymetry_list), str(input_bathymetry_list) + ".yxz"),
            header=None,
            names=["y", "x", "depth"], sep=" ")

        depth_array.loc[depth_array["depth"] == no_data_value, "depth"] = np.nan

        depth_array_min = depth_array['depth'].min()
        depth_array_max = depth_array['depth'].max()
        depth_array_x_min = depth_array['x'].min()
        depth_array_x_max = depth_array['x'].max()
        depth_array_y_min = depth_array['y'].min()
        depth_array_y_max = depth_array['y'].max()

        tri_value_list = []

        if depth_array_min < 0:
            depth_array["depth"] *= -1
            depth_array_min = depth_array['depth'].min()
            depth_array_max = depth_array['depth'].max()

        def build_env_array(location, bname, z_min, z_max, y_min, y_max, x_min, x_max):
            # 1 Extract list of files.
            env_files = []
            env_values = []
            env.workspace = location
            rasterlist1 = arcpy.ListRasters("*")
            for f in rasterlist1:
                f2 = float(f.replace(bname, ""))
                env_values.append(f2)
                env_files.append(f)
            else:
                pass

            env_file_list = sorted(zip(env_files, env_values), key=lambda tup: tup[1])

            # 2 Get xyz values needed to build array
            xy_coords = pd.read_pickle(os.path.join(location, "xy_coords.yxz"))

            y_min = y_min - (input_environment0_cs)
            y_max = y_max + (input_environment0_cs)
            x_min = x_min - (input_environment0_cs)
            x_max = x_max + (input_environment0_cs)

            x_vals = np.unique(xy_coords["x"])
            y_vals = np.unique(xy_coords["y"])

            if input_environment0_cs_x_min > x_min and input_environment0_cs_y_min > y_min:
                x_min = input_environment0_cs_x_min
                y_min = input_environment0_cs_y_min
            elif input_environment0_cs_x_min > x_min:
                x_min = input_environment0_cs_x_min
            elif input_environment0_cs_y_min > y_min:
                y_min = input_environment0_cs_y_min

            if input_environment0_cs_x_max < x_max and input_environment0_cs_y_max < y_max:
                x_max = input_environment0_cs_x_max
                y_max = input_environment0_cs_y_max
            elif input_environment0_cs_x_max < x_max:
                x_max = input_environment0_cs_x_max
            elif input_environment0_cs_y_max < y_max:
                y_max = input_environment0_cs_y_max
                if verbose_mode:
                    print("4")

            if verbose_mode:
                print("y min: " + str(y_min))
                print("y max: " + str(y_max))

            x_vals = [i for i in x_vals if i > x_min and i < x_max]
            y_vals = [i for i in y_vals if i > y_min and i < y_max]

            # Some times we may have bathymetric points outside of the environmental
            # layer, and as such, we need to work around this. It only works for the
            # situation where y_min > input_environment0_cs_y_max, as the SplitRaster code
            # will not usually have a problem at the input_environment0_cs_y_min region.
            if not y_vals:
                y_min = input_environment0_cs_y_max - (input_environment0_cs * 3)
                y_vals = [i for i in np.unique(xy_coords["y"]) if i > y_min and i < y_max]

            x_vals_min = min(x_vals)
            x_vals_max = max(x_vals)
            y_vals_min = min(y_vals)
            y_vals_max = max(y_vals)

            if y_vals_min == y_vals_max:
                y_min = y_vals_min - (input_environment0_cs * 3)
                y_vals = [i for i in np.unique(xy_coords["y"]) if i > y_min and i < y_max]
                y_vals_min = min(y_vals)
                y_vals_max = max(y_vals)

            temp_name = []
            temp_depth = []

            counter = 1
            key = -1

            for t_pair in env_file_list:
                key += 1
                name, depth = t_pair
                if z_min <= depth <= z_max:
                    key2 = key
                    if counter == 1:
                        v1, v2 = env_file_list[key - 1]
                        temp_name.append(v1)
                        temp_depth.append(v2)
                        del v1, v2
                    temp_name.append(name)
                    temp_depth.append(depth)
                    counter += 1

            if z_min == z_max:
                key = -1
                for t_pair in env_file_list:
                    key += 1
                    name, depth = t_pair
                    if depth <= z_min:
                        key2 = key
                        if counter == 1:
                            v1, v2 = env_file_list[key]
                            temp_name.append(v1)
                            temp_depth.append(v2)
                            del v1, v2
                        counter += 1
            if "key2" in locals():
                if key2 < key:
                    v1, v2 = env_file_list[key2 + 1]
                    temp_name.append(v1)
                    temp_depth.append(v2)
                    del v1, v2
            else:
                key = -1
                for t_pair in env_file_list:
                    key += 1
                    name, depth = t_pair
                    if depth >= z_min:
                        if counter == 1:
                            key2 = key
                            v1, v2 = env_file_list[key - 1]
                            temp_name.append(v1)
                            temp_depth.append(v2)
                            del v1, v2
                        counter += 1

                if key2 < key:
                    v1, v2 = env_file_list[key2 + 1]
                    temp_name.append(v1)
                    temp_depth.append(v2)
                    del v1, v2

            temp_depth_reverse = temp_depth[::-1]

            z_vals = np.unique(np.asarray(temp_depth_reverse, dtype=float))

            if verbose_mode:
                print "zval len=" + str(z_vals.size) + " , z_vals=" + str(z_vals)

            def pad(data):
                bad_indexes = np.isnan(data)
                good_indexes = np.logical_not(bad_indexes)
                good_data = data[good_indexes]
                interpolated = np.interp(bad_indexes.nonzero()[0], good_indexes.nonzero()[0], good_data)
                data[bad_indexes] = interpolated
                return data

            data = np.array([np.flipud(pad(arcpy.RasterToNumPyArray(os.path.join(location, bname + "%f" % f),
                                                                    arcpy.Point(x_min, y_min), len(x_vals),
                                                                    len(y_vals), nodata_to_value=np.nan))) for f in
                             temp_depth])
            data = data.T

            file = open(os.path.join(output_directory, "Outputs",
                                     "env_" + str(input_bathymetry_split_i) + ".csv"), 'w')
            for x in x_vals:
                for y in y_vals:
                    file.write(str(x) + ", " + str(y) + ", " + '\n')

            file.close()
            if verbose_mode:
                print "1: x- " + str(x_vals) + " y- " + str(y_vals) + " z- " + str(z_vals)
                print "here"
                print data.size
                print len(x_vals)
                print len(y_vals)
                print len(z_vals)

            rgi = RegularGridInterpolator((x_vals, y_vals, z_vals), data, method='linear')
            return rgi, y_vals_min, y_vals_max, x_vals_min, x_vals_max

        # Build extents for subsetting the environemtnal grid

        if depth_array_min == depth_array_min_comp and depth_array_max == depth_array_max_comp:
            arcpy.AddMessage(
                "Skipped building environment value array for " + str(input_bathymetry_split_i))
        else:
            depth_array_min_comp <= depth_array_min
            depth_array_max_comp >= depth_array_max
            arcpy.AddMessage("Building environment value array: " + str(input_bathymetry_split_i))

            # deal with single row/column data
            if depth_array_y_min == depth_array_y_max:
                depth_array_y_min = depth_array_y_min - (input_environment0_cs * 3)
                depth_array_y_max = depth_array_y_max + (input_environment0_cs * 3)

            if depth_array_x_min == depth_array_x_max:
                depth_array_x_min = depth_array_x_min - (input_environment0_cs * 3)
                depth_array_x_max = depth_array_x_max + (input_environment0_cs * 3)

            if verbose_mode:
                print "2: env- " + str(input_environment) + " envname- " + str(input_environment_name) + " depthmin- " + \
                  str(depth_array_min) + " depthmax- " + str(depth_array_max) + " depth_array_y_min- " + str(depth_array_y_min) \
                  + " depth_array_y_max- " + str(depth_array_y_max) + " depth_array_x_min- " + str(depth_array_x_min) \
                  + " depth_array_x_max- " + str(depth_array_x_max)

            rgi, y_vals_min, y_vals_max, x_vals_min, x_vals_max = build_env_array(input_environment,
                                                                                  input_environment_name,
                                                                                  depth_array_min,
                                                                                  depth_array_max,
                                                                                  depth_array_y_min,
                                                                                  depth_array_y_max,
                                                                                  depth_array_x_min,
                                                                                  depth_array_x_max)

        depth_array = depth_array.iloc[::-1]
        arcpy.AddMessage("Trilinearly interpolating: " + str(input_bathymetry_split_i))

        # Check for edge issue
        if depth_array_x_min < x_vals_min or depth_array_x_max > x_vals_max or depth_array_y_min < y_vals_min or depth_array_y_max > y_vals_max:
            edge = True
            if verbose_mode:
                arcpy.AddMessage(str(input_bathymetry_split_i) + " is edge")
        else:
            edge = False
            if verbose_mode:
                arcpy.AddMessage(str(input_bathymetry_split_i) + " is not edge")

        f = open(os.path.join(output_directory, "Outputs",
                              "tri_" + str(input_bathymetry_split_i) + ".csv"), 'w')

        if not edge:
            for index, row in depth_array.iterrows():
                if verbose_mode:
                    print str(index) + " " + str(row)

                if not np.isnan(row["depth"]) and row["depth"] >= 5500.:
                    tri_value_list.append(rgi((row["x"], row["y"], 5499.0)))
                    f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                        rgi((row["x"], row["y"], 5499.0))) + ', Flag 4\n')
                elif not np.isnan(row["depth"]) and row["depth"] < 0.:
                    tri_value_list.append(rgi((row["x"], row["y"], 0.0)))
                    f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                        rgi((row["x"], row["y"], 0.0))) + ', Flag 4\n')
                elif np.isnan(row["depth"]):
                    tri_value_list.append("-9999")
                    f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + "-9999" + ", Flag 3\n")
                elif not np.isnan(row["depth"]):
                    tri_value_list.append(rgi((row["x"], row["y"], row["depth"])))
                    f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                        rgi((row["x"], row["y"], row["depth"]))) + ", Flag 2" + '\n')
        elif edge:
            # Deal with edge effect
            for index, row in depth_array.iterrows():
                if not np.isnan(row["depth"]) and row["depth"] >= 5500.:
                    # In appropriate space
                    if row["x"] > x_vals_min and row["x"] < x_vals_max and row["y"] > y_vals_min and row[
                        "y"] < y_vals_max:
                        tri_value_list.append(rgi((row["x"], row["y"], 5499.0)))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((row["x"], row["y"], 5499.0))) + ", OKDeep" + '\n')
                    # Bottom Left
                    elif row["x"] < x_vals_min and row["y"] < y_vals_min:
                        tri_value_list.append(rgi((x_vals_min, y_vals_min, 5499.0)))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((x_vals_min, y_vals_min, 5499.0))) + ", BotLeftDeep" + '\n')
                    # Upper Left
                    elif row["x"] < x_vals_min and row["y"] > y_vals_max:
                        tri_value_list.append(rgi((x_vals_min, y_vals_max, 5499.0)))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((row["x"], row["y"], 5499.0))) + ", UpLeftDeep" + '\n')
                    # Bottom Right
                    elif row["x"] > x_vals_max and row["y"] < y_vals_min:
                        tri_value_list.append(rgi((x_vals_max, y_vals_min, 5499.0)))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((x_vals_max, y_vals_min, 5499.0))) + ", BotRightDeep" + '\n')
                    # Upper Right
                    elif row["x"] > x_vals_max and row["y"] > y_vals_max:
                        tri_value_list.append(rgi((x_vals_max, y_vals_max, 5499.0)))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((x_vals_max, y_vals_max, 5499.0))) + ", UpRightDeep" + '\n')
                    # Top
                    elif row["y"] > y_vals_max:
                        tri_value_list.append(rgi((row["x"], y_vals_max, 5499.0)))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((row["x"], y_vals_max, 5499.0))) + ", TopDeep" + '\n')
                    # Bottom
                    elif row["y"] < y_vals_min:
                        tri_value_list.append(rgi((row["x"], y_vals_min, 5499.0)))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((row["x"], y_vals_min, 5499.0))) + ", BotDeep" + '\n')
                    # Left
                    elif row["x"] < x_vals_min:
                        tri_value_list.append(rgi((x_vals_min, row["y"], 5499.0)))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((x_vals_min, row["y"], 5499.0))) + ", LeftDeep" + '\n')
                    # Right
                    elif row["x"] > x_vals_max:
                        tri_value_list.append(rgi((x_vals_max, row["y"], 5499.0)))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((x_vals_max, row["y"], 5499.0))) + ", RightDeep" + '\n')
                elif np.isnan(row["depth"]):
                    tri_value_list.append("-9999")
                    f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + "-9999" + ", Flagged 1" + '\n')
                elif not np.isnan(row["depth"]):
                    # In appropriate space
                    if row["x"] > x_vals_min and row["x"] < x_vals_max and row["y"] > y_vals_min and row[
                        "y"] < y_vals_max:
                        tri_value_list.append(rgi((row["x"], row["y"], row["depth"])))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((row["x"], row["y"], row["depth"]))) + ", OK, " + str(row["depth"]) + '\n')
                    # Bottom Left
                    elif row["x"] < x_vals_min and row["y"] < y_vals_min:
                        tri_value_list.append(rgi((x_vals_min, y_vals_min, row["depth"])))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((x_vals_min, y_vals_min, row["depth"]))) + ", BotLeft" + '\n')
                    # Upper Left
                    elif row["x"] < x_vals_min and row["y"] > y_vals_max:
                        tri_value_list.append(rgi((x_vals_min, y_vals_max, row["depth"])))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((x_vals_min, y_vals_max, row["depth"]))) + ", UpLeft" + '\n')
                    # Bottom Right
                    elif row["x"] > x_vals_max and row["y"] < y_vals_min:
                        tri_value_list.append(rgi((x_vals_max, y_vals_min, row["depth"])))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((x_vals_max, y_vals_min, row["depth"]))) + ", BotRight" + '\n')
                    # Upper Right
                    elif row["x"] > x_vals_max and row["y"] > y_vals_max:
                        tri_value_list.append(rgi((x_vals_max, y_vals_max, row["depth"])))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((x_vals_max, y_vals_max, row["depth"]))) + ", UpRight" + '\n')
                    # Top
                    elif row["y"] > y_vals_max:
                        tri_value_list.append(rgi((row["x"], y_vals_max, row["depth"])))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((row["x"], y_vals_max, row["depth"]))) + ", Top" + '\n')
                    # Bottom
                    elif row["y"] < y_vals_min:
                        tri_value_list.append(rgi((row["x"], y_vals_min, row["depth"])))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((row["x"], y_vals_min, row["depth"]))) + ", Bottom" + '\n')
                    # Left
                    elif row["x"] < x_vals_min:
                        tri_value_list.append(rgi((x_vals_min, row["y"], row["depth"])))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((x_vals_min, row["y"], row["depth"]))) + ", Left" + '\n')
                    # Right
                    elif row["x"] > x_vals_max:
                        tri_value_list.append(rgi((x_vals_max, row["y"], row["depth"])))
                        f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            rgi((x_vals_max, row["y"], row["depth"]))) + ", Right" + '\n')
                else:
                    f.write(str(row["x"]) + ", " + str(row["y"]) + ", 100000" + ', Flag 5\n')

        f.close()

        # Get headers
        with open(os.path.join(output_directory,
                               str(input_bathymetry_split_i) + ".asc")) as myfile:
            header_ascii = myfile.readlines()[0:6]  # put here the interval you want
        # arcpy.AddMessage("ASCII Header: " + str(header_ascii))

        # Write headers to new file
        with open(os.path.join(output_directory, "Outputs",
                               str(input_bathymetry_split_i) + ".asc"), 'a') as ascii_file:
            for wr in header_ascii:
                ascii_file.write(wr)

        with open(os.path.join(output_directory, "Outputs",
                               str(input_bathymetry_split_i) + ".asc"), 'a') as ascii_file:
            for tri in tri_value_list:
                ascii_file.write(str(tri) + " ")

        arcpy.AddMessage("Completed block: " + str(input_bathymetry_split_i))

        del tri_value_list, tri, depth_array, depth_array_min, depth_array_max, input_bathymetry_split
    else:
        arcpy.AddMessage("Skipping empty or previously completed raster: " + str(input_bathymetry_split_i))
        del input_bathymetry_split

    return


# Stage 3 - Mosaic the chunks back into the extent of the original bathymetry

def mpchunk(output_directory, input_bathymetry_cs, output_list_chunk):
    arcpy.AddMessage("Processing chunk: " + str(output_list_chunk[0]))
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

    arcpy.AddMessage("Completed chunk: " + str(output_list_chunk[0]))

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


if __name__ == '__main__':

    ## Temperature
    input_environment = r"D:\sponges\world-ocean-atlas\Inputs\temperature\gridded_0-2\Projected"
    environment_name = r"t_an"
    cpu_cores_used = "44"
    test_mode = False
    verbose_mode = False

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\etopo2oceanp_sr"
    output_directory = r"D:\sponges\bathymetry\Global\etopo2oceanp_temp"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\gebco08p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\gebco08p_temp"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\gebco14p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\gebco14p_temp"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\srtm30p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\srtm30p_temp"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\srtm15p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\srtm15p_temp"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    ## Salinity
    input_environment = r"D:\sponges\world-ocean-atlas\Inputs\salinity\Projected"
    environment_name = r"s_an"
    cpu_cores_used = "44"
    test_mode = False
    verbose_mode = False

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\etopo2oceanp_sr"
    output_directory = r"D:\sponges\bathymetry\Global\etopo2oceanp_sal"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\gebco08p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\gebco08p_sal"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\gebco14p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\gebco14p_sal"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)
    
    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\srtm30p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\srtm30p_sal"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\srtm15p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\srtm15p_sal"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    ## Apparent oxygen utilisation
    input_environment = r"D:\sponges\world-ocean-atlas\Inputs\aoxu\gridded_0-5_bat\Projected"
    environment_name = r"a_an"
    cpu_cores_used = "44"
    test_mode = False
    verbose_mode = False

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\etopo2oceanp_sr"
    output_directory = r"D:\sponges\bathymetry\Global\etopo2oceanp_aoxu"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode, verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\gebco08p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\gebco08p_aoxu"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode, verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\gebco14p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\gebco14p_aoxu"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode, verbose_mode)
    
    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\srtm30p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\srtm30p_aoxu"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\srtm15p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\srtm15p_aoxu"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    ## Dissolved oxygen
    input_environment = r"D:\sponges\world-ocean-atlas\Inputs\dissolvedO2\extracted\Projected"
    environment_name = r"o_an"
    cpu_cores_used = "44"
    test_mode = False
    verbose_mode = False

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\etopo2oceanp_sr"
    output_directory = r"D:\sponges\bathymetry\Global\etopo2oceanp_diso2"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode, verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\gebco08p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\gebco08p_diso2"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode, verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\gebco14p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\gebco14p_diso2"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode, verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\srtm30p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\srtm30p_diso2"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\srtm15p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\srtm15p_diso2"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)
    
    ## Nitrate
    input_environment = r"D:\sponges\world-ocean-atlas\Inputs\nitrate\Projected"
    environment_name = r"n_an"
    cpu_cores_used = "44"
    test_mode = False
    verbose_mode = False

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\etopo2oceanp_sr"
    output_directory = r"D:\sponges\bathymetry\Global\etopo2oceanp_nit"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\gebco08p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\gebco08p_nit"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\gebco14p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\gebco14p_nit"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)
    
    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\srtm30p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\srtm30p_nit"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\srtm15p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\srtm15p_nit"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    ## Phosphate
    input_environment = r"D:\sponges\world-ocean-atlas\Inputs\phosphate\Projected"
    environment_name = r"p_an"
    cpu_cores_used = "44"
    test_mode = False
    verbose_mode = False

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\etopo2oceanp_sr"
    output_directory = r"D:\sponges\bathymetry\Global\etopo2oceanp_phos"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\gebco08p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\gebco08p_phos"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\gebco14p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\gebco14p_phos"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\srtm30p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\srtm30p_phos"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\srtm15p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\srtm15p_phos"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    ## Silicate
    input_environment = r"D:\sponges\world-ocean-atlas\Inputs\silicate\Projected"
    environment_name = r"i_an"
    cpu_cores_used = "44"
    test_mode = False
    verbose_mode = False

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\etopo2oceanp_sr"
    output_directory = r"D:\sponges\bathymetry\Global\etopo2oceanp_sil"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\gebco08p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\gebco08p_sil"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\gebco14p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\gebco14p_sil"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)
    
    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\srtm30p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\srtm30p_sil"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)

    input_bathymetry_folder = r"D:\sponges\bathymetry\Global\srtm15p_sr"
    output_directory = r"D:\sponges\bathymetry\Global\srtm15p_sil"

    globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, test_mode,
            verbose_mode)
