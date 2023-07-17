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

import sys, os, arcpy, time, multiprocessing, csv, gc
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import time
from arcpy import env
from arcpy.sa import *
import pandas as pd
from functools import partial
import subprocess


def globenv(input_bathymetry_folder, input_environment, environment_name, output_directory, cpu_cores_used, chunk_mode, test_mode, verbose_mode):

    verbose_mode = verbose_mode
    environment_name = environment_name

    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("Spatial")

    print "GlobENV - Trilinear interpolation of environmental data onto bathymetry"
    print "... Processing " + str(input_bathymetry_folder)

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

    print "... There are " + str(len(input_environment_depth)) + " environmental layers"

    if len(input_environment_depth) < 2:
        print "Error: Too few environmental layers for Trilinear interpolation"
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
        print "... Both depth and environmental data are in valid projection"
        if input_bathymetry_sr.name == input_environment0_sr.name:
            pass
        else:
            print "Error: Depth and environmental data are in different projections, must be identical"
            sys.exit(0)
        pass
    else:
        print "Error: Both depth and environmental data need to be in projected coordinate system"
        sys.exit(0)

    if not os.path.exists(os.path.join(output_directory, "Outputs")):
        os.makedirs(os.path.join(output_directory, "Outputs"))

    env.workspace = input_bathymetry_folder
    input_bathymetry_list = arcpy.ListRasters("sp*", "GRID")

    print "... There are " + str(len(input_bathymetry_list)) + " depth blocks to process."

    input_bathymetry_list_done = []

    for i in input_bathymetry_list:
        try:
            r = arcpy.Raster(os.path.join(input_bathymetry_folder, str(i)))
        except:
            r = False
            print "... Issue reading raster " + str(i)

        del r
        if not os.path.exists(os.path.join(output_directory, "Outputs", str(i) + ".asc")):
            input_bathymetry_list_done.append(i)

    if len(input_bathymetry_list_done) > 0:
        if test_mode:
            print "... Test mode: Processing in a single thread"
            for i in input_bathymetry_list_done:
                mpprocess(output_directory, input_environment_depth, input_environment0_cs, input_environment,
                           input_environment_name, input_environment0_cs_x_min, input_environment0_cs_y_min,
                           input_environment0_cs_x_max, input_environment0_cs_y_max, input_bathymetry_folder, verbose_mode, i)
        else:
            print "... Parallel mode: There are " + str(len(input_bathymetry_list_done)) + " depth left to process."

            if int(len(input_bathymetry_list_done)) < int(cpu_cores_used):
                pool = multiprocessing.Pool(int(len(input_bathymetry_list_done)))
                print "... Will use " + str(int(len(input_bathymetry_list_done))) + " cores for processing"
            else:
                pool = multiprocessing.Pool(int(cpu_cores_used))
                print "... Will use " + str(cpu_cores_used) + " cores for processing"

            func = 1
            ans = pool.map(partial(mpprocess, output_directory, input_environment_depth, input_environment0_cs, input_environment,
                           input_environment_name, input_environment0_cs_x_min, input_environment0_cs_y_min,
                           input_environment0_cs_x_max,
                           input_environment0_cs_y_max, input_bathymetry_folder, verbose_mode), input_bathymetry_list_done)
            pool.close()
            pool.join()

    else:
        print "... All rasters processed"


    if chunk_mode == 'gdal':
        env.workspace = os.path.join(output_directory, "Outputs")
        output_list = arcpy.ListRasters("*", "ALL")
        gdalchunk(output_directory, output_list)

    elif chunk_mode == 'arcpy':

        print "... Mosaicking using arcpy, very slow"

        env.workspace = os.path.join(output_directory, "Outputs")
        output_list = arcpy.ListRasters("*", "ALL")

        if not os.path.exists(os.path.join(output_directory, "Outputs_C")):
            os.makedirs(os.path.join(output_directory, "Outputs_C"))

        if len(output_list) > 200:
            output_list_chunk = [output_list[x:x + 20] for x in xrange(0, len(output_list), 20)]
        else:
            output_list_chunk = [output_list[x:x + 5] for x in xrange(0, len(output_list), 5)]

        print "... Building pool to mosaic " + str(len(output_list)) + " rasters, " + "in " + str(
            len(output_list_chunk)) + " chunks."

        if int(len(output_list_chunk)) < int(cpu_cores_used):
            pool2 = multiprocessing.Pool(int(len(output_list_chunk)))
            print "... Will use " + str(int(len(output_list_chunk))) + " cores for mosaic processing"
        else:
            pool2 = multiprocessing.Pool(int(cpu_cores_used))
            print "... Will use " + str(cpu_cores_used) + " cores for mosaic processing"

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

        print "... Mosaicking " + str(len(chunk_list_a)) + " chunks together."

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

        print "... GlobENV complete!"

def mpprocess(output_directory, input_environment_depth, input_environment0_cs, input_environment,
                       input_environment_name, input_environment0_cs_x_min, input_environment0_cs_y_min,
                       input_environment0_cs_x_max,
                       input_environment0_cs_y_max, input_bathymetry_folder, verbose_mode, input_bathymetry_list):
    try:
        arcpy.CheckOutExtension("Spatial")
        arcpy.env.overwriteOutput = True
        depth_array_min_comp = -9999
        depth_array_max_comp = -9999

        if not os.path.exists(os.path.join(output_directory, "1_Temp", str(input_bathymetry_list))):
            os.makedirs(os.path.join(output_directory, "1_Temp", str(input_bathymetry_list)))

        env.workspace = os.path.join(output_directory, "1_Temp", str(input_bathymetry_list))
        env.scratchWorkspace = os.path.join(output_directory, "1_Temp", str(input_bathymetry_list))

        log_file("GlobENV Logfile\n\n", output_directory, input_bathymetry_list, verbose_mode)
        log_file("Bathymetry: " + str(input_bathymetry_list) + "\n", output_directory, input_bathymetry_list, verbose_mode)
        log_file("Environment name: " + str(input_environment_name) + "\n", output_directory, input_bathymetry_list, verbose_mode)

        try:
            arcpy.RasterToASCII_conversion(os.path.join(input_bathymetry_folder, input_bathymetry_list),
                                           os.path.join(output_directory,
                                                        str(
                                                            input_bathymetry_list) + ".asc"))
            input_bathymetry_split = arcpy.Raster(os.path.join(output_directory,
                                                               str(
                                                                   input_bathymetry_list) + ".asc"))
            input_bathymetry_split_i = input_bathymetry_list
        except:
            log_file.write("Unable to read bathymetric layer: " + str(input_bathymetry_list))
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
                xy_coords = pd.read_pickle(os.path.join(location, "xy_coords.pkl"))

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

                log_file("y min: " + str(y_min) + "\n", output_directory, input_bathymetry_list, verbose_mode)
                log_file("y max: " + str(y_max) + "\n", output_directory, input_bathymetry_list, verbose_mode)
                log_file("x min: " + str(x_min) + "\n", output_directory, input_bathymetry_list, verbose_mode)
                log_file("x max: " + str(x_max) + "\n", output_directory, input_bathymetry_list, verbose_mode)
                log_file("z min: " + str(z_min) + "\n", output_directory, input_bathymetry_list, verbose_mode)
                log_file("z max: " + str(z_max) + "\n", output_directory, input_bathymetry_list, verbose_mode)

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

                env_name = []
                env_depth = []

                for pair in env_file_list:
                    env_name_i, env_depth_i = pair
                    env_name.append(env_name_i)
                    env_depth.append(env_depth_i)

                log_file("Env array (items): " + str(env_depth) + "\n", output_directory, input_bathymetry_list, verbose_mode)
                log_file("Env array (length): " + str(len(env_depth)) + "\n", output_directory, input_bathymetry_list, verbose_mode)

                # Get Indices of location in list
                env_min_array_depth = min(range(len(env_depth)), key=lambda i: abs(env_depth[i] - z_min))
                env_max_array_depth = min(range(len(env_depth)), key=lambda i: abs(env_depth[i] - z_max))

                log_file("env array zmin (pre pad): " + str(env_depth[env_min_array_depth]) + "\n", output_directory, input_bathymetry_list, verbose_mode)
                log_file("env array zmax (pre pad): " + str(env_depth[env_max_array_depth]) + "\n", output_directory, input_bathymetry_list, verbose_mode)

                if env_min_array_depth > 0:
                    env_min_array_depth = env_min_array_depth - 1
                    log_file("env array zmin (post pad): " + str(env_depth[env_min_array_depth]) + "\n", output_directory,
                             input_bathymetry_list, verbose_mode)

                if env_max_array_depth < len(env_depth) - 1:
                    env_max_array_depth = env_max_array_depth + 2
                    log_file("Pad value: +2" + "\n",
                             output_directory,
                             input_bathymetry_list, verbose_mode)
                    env_array_depth = env_depth[env_min_array_depth:env_max_array_depth]


                if env_max_array_depth >= len(env_depth) - 1:
                    env_max_array_depth = len(env_depth) - 1
                    log_file("Pad value: -1" + "\n",
                             output_directory,
                             input_bathymetry_list, verbose_mode)
                    env_array_depth = env_depth[env_min_array_depth:]

                log_file("env array zmax (index value): " + str(env_max_array_depth) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)

                log_file("env array zmax (post pad): " + str(env_depth[env_max_array_depth]) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)

                log_file("Environment depth array (values list): " + str(env_array_depth) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)
                env_array_depth_reverse = env_array_depth[::-1]
                log_file("Environment depth array (values list reverse): " + str(env_array_depth) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)
                z_vals = np.unique(np.asarray(env_array_depth_reverse, dtype=float))

                log_file("Environment depth array (number): " + str(z_vals.size) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)

                log_file("Environment depth array (values): " + str(z_vals) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)

                def pad(data):
                    bad_indexes = np.isnan(data)
                    good_indexes = np.logical_not(bad_indexes)
                    good_data = data[good_indexes]
                    interpolated = np.interp(bad_indexes.nonzero()[0], good_indexes.nonzero()[0], good_data)
                    data[bad_indexes] = interpolated
                    return data

                data = np.array([np.flipud(pad(arcpy.RasterToNumPyArray(os.path.join(location, bname + "%f" % f),
                                                                        arcpy.Point(x_min, y_min), len(x_vals),
                                                                        len(y_vals), nodata_to_value=np.nan))) for f in env_array_depth_reverse])
                data = data.T

                file = open(os.path.join(output_directory, "Outputs",
                                         "env_" + str(input_bathymetry_split_i) + ".csv"), 'w')
                for x in x_vals:
                    for y in y_vals:
                        file.write(str(x) + ", " + str(y) + ", " + '\n')

                file.close()

                log_file("Environment TRI array (np shape): " + str(data.shape) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)

                rgi = RegularGridInterpolator((x_vals, y_vals, z_vals), data, method='linear')
                return rgi, y_vals_min, y_vals_max, x_vals_min, x_vals_max

            # Build extents for subsetting the environemtnal grid

            if depth_array_min == depth_array_min_comp and depth_array_max == depth_array_max_comp:
                log_file("Skipped building environment value array: " + str(input_bathymetry_split_i) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)
            else:
                depth_array_min_comp <= depth_array_min
                depth_array_max_comp >= depth_array_max

                log_file("Building environment value array: " + str(input_bathymetry_split_i) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)

                # deal with single row/column data
                if depth_array_y_min == depth_array_y_max:
                    depth_array_y_min = depth_array_y_min - (input_environment0_cs * 3)
                    depth_array_y_max = depth_array_y_max + (input_environment0_cs * 3)

                if depth_array_x_min == depth_array_x_max:
                    depth_array_x_min = depth_array_x_min - (input_environment0_cs * 3)
                    depth_array_x_max = depth_array_x_max + (input_environment0_cs * 3)

                log_file("\nTRI Parameters" + "\n", output_directory, input_bathymetry_list, verbose_mode)
                log_file("Environment data: " + str(input_environment) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)
                log_file("Environment name: " + str(input_environment_name) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)
                log_file("Depth min: " + str(depth_array_min) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)
                log_file("Depth max: " + str(depth_array_max) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)
                log_file("Depth array y min: " + str(depth_array_y_min) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)
                log_file("Depth array y max: " + str(depth_array_y_max) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)
                log_file("Depth array xmin: " + str(depth_array_x_min) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)
                log_file("Depth array x max: " + str(depth_array_max) + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)

                rgi, y_vals_min, y_vals_max, x_vals_min, x_vals_max = build_env_array(input_environment,
                                                                                      input_environment_name,
                                                                                      depth_array_min,
                                                                                      depth_array_max,
                                                                                      depth_array_y_min,
                                                                                      depth_array_y_max,
                                                                                      depth_array_x_min,
                                                                                      depth_array_x_max)

            depth_array = depth_array.iloc[::-1]

            log_file("Trilinearly interpolating: " + str(input_bathymetry_split_i) + "\n", output_directory,
                     input_bathymetry_list, verbose_mode)

            # Check for edge issues
            if depth_array_x_min < x_vals_min or depth_array_x_max > x_vals_max or depth_array_y_min < y_vals_min or depth_array_y_max > y_vals_max:
                edge = True
                log_file("Edge status: " +str(input_bathymetry_split_i) + " is edge" + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)
            else:
                edge = False
                log_file("Edge status: " +str(input_bathymetry_split_i) + " is not edge" + "\n", output_directory,
                         input_bathymetry_list, verbose_mode)

            f = open(os.path.join(output_directory, "Outputs",
                                  "tri_" + str(input_bathymetry_split_i) + ".csv"), 'w')

            if not edge:
                for index, row in depth_array.iterrows():
                    # if verbose_mode:
                    #     print str(index) + " " + str(row)

                    if not np.isnan(row["depth"]) and row["depth"] >= 5500.:
                        tri_value_list.append(rgi((row["x"], row["y"], 5499.0)))
                        # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                        #     rgi((row["x"], row["y"], 5499.0))) + ', Flag 4\n')
                    elif not np.isnan(row["depth"]) and row["depth"] < 0.:
                        tri_value_list.append(rgi((row["x"], row["y"], 0.0)))
                        # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                        #     rgi((row["x"], row["y"], 0.0))) + ', Flag 4\n')
                    elif np.isnan(row["depth"]):
                        tri_value_list.append("-9999")
                        # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + "-9999" + ", Flag 3\n")
                    elif not np.isnan(row["depth"]):
                        tri_value_list.append(rgi((row["x"], row["y"], row["depth"])))
                        # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                        #     rgi((row["x"], row["y"], row["depth"]))) + ", Flag 2" + '\n')
            elif edge:
                # Deal with edge effect
                for index, row in depth_array.iterrows():
                    if not np.isnan(row["depth"]) and row["depth"] >= 5500.:
                        # In appropriate space
                        if row["x"] > x_vals_min and row["x"] < x_vals_max and row["y"] > y_vals_min and row[
                            "y"] < y_vals_max:
                            tri_value_list.append(rgi((row["x"], row["y"], 5499.0)))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((row["x"], row["y"], 5499.0))) + ", OKDeep" + '\n')
                        # Bottom Left
                        elif row["x"] < x_vals_min and row["y"] < y_vals_min:
                            tri_value_list.append(rgi((x_vals_min, y_vals_min, 5499.0)))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((x_vals_min, y_vals_min, 5499.0))) + ", BotLeftDeep" + '\n')
                        # Upper Left
                        elif row["x"] < x_vals_min and row["y"] > y_vals_max:
                            tri_value_list.append(rgi((x_vals_min, y_vals_max, 5499.0)))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((row["x"], row["y"], 5499.0))) + ", UpLeftDeep" + '\n')
                        # Bottom Right
                        elif row["x"] > x_vals_max and row["y"] < y_vals_min:
                            tri_value_list.append(rgi((x_vals_max, y_vals_min, 5499.0)))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((x_vals_max, y_vals_min, 5499.0))) + ", BotRightDeep" + '\n')
                        # Upper Right
                        elif row["x"] > x_vals_max and row["y"] > y_vals_max:
                            tri_value_list.append(rgi((x_vals_max, y_vals_max, 5499.0)))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((x_vals_max, y_vals_max, 5499.0))) + ", UpRightDeep" + '\n')
                        # Top
                        elif row["y"] > y_vals_max:
                            tri_value_list.append(rgi((row["x"], y_vals_max, 5499.0)))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((row["x"], y_vals_max, 5499.0))) + ", TopDeep" + '\n')
                        # Bottom
                        elif row["y"] < y_vals_min:
                            tri_value_list.append(rgi((row["x"], y_vals_min, 5499.0)))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((row["x"], y_vals_min, 5499.0))) + ", BotDeep" + '\n')
                        # Left
                        elif row["x"] < x_vals_min:
                            tri_value_list.append(rgi((x_vals_min, row["y"], 5499.0)))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((x_vals_min, row["y"], 5499.0))) + ", LeftDeep" + '\n')
                        # Right
                        elif row["x"] > x_vals_max:
                            tri_value_list.append(rgi((x_vals_max, row["y"], 5499.0)))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((x_vals_max, row["y"], 5499.0))) + ", RightDeep" + '\n')
                    elif np.isnan(row["depth"]):
                        tri_value_list.append("-9999")
                        # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + "-9999" + ", Flagged 1" + '\n')
                    elif not np.isnan(row["depth"]):
                        # In appropriate space
                        if row["x"] > x_vals_min and row["x"] < x_vals_max and row["y"] > y_vals_min and row[
                            "y"] < y_vals_max:
                            tri_value_list.append(rgi((row["x"], row["y"], row["depth"])))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((row["x"], row["y"], row["depth"]))) + ", OK, " + str(row["depth"]) + '\n')
                        # Bottom Left
                        elif row["x"] < x_vals_min and row["y"] < y_vals_min:
                            tri_value_list.append(rgi((x_vals_min, y_vals_min, row["depth"])))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((x_vals_min, y_vals_min, row["depth"]))) + ", BotLeft" + '\n')
                        # Upper Left
                        elif row["x"] < x_vals_min and row["y"] > y_vals_max:
                            tri_value_list.append(rgi((x_vals_min, y_vals_max, row["depth"])))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((x_vals_min, y_vals_max, row["depth"]))) + ", UpLeft" + '\n')
                        # Bottom Right
                        elif row["x"] > x_vals_max and row["y"] < y_vals_min:
                            tri_value_list.append(rgi((x_vals_max, y_vals_min, row["depth"])))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((x_vals_max, y_vals_min, row["depth"]))) + ", BotRight" + '\n')
                        # Upper Right
                        elif row["x"] > x_vals_max and row["y"] > y_vals_max:
                            tri_value_list.append(rgi((x_vals_max, y_vals_max, row["depth"])))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((x_vals_max, y_vals_max, row["depth"]))) + ", UpRight" + '\n')
                        # Top
                        elif row["y"] > y_vals_max:
                            tri_value_list.append(rgi((row["x"], y_vals_max, row["depth"])))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((row["x"], y_vals_max, row["depth"]))) + ", Top" + '\n')
                        # Bottom
                        elif row["y"] < y_vals_min:
                            tri_value_list.append(rgi((row["x"], y_vals_min, row["depth"])))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((row["x"], y_vals_min, row["depth"]))) + ", Bottom" + '\n')
                        # Left
                        elif row["x"] < x_vals_min:
                            tri_value_list.append(rgi((x_vals_min, row["y"], row["depth"])))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((x_vals_min, row["y"], row["depth"]))) + ", Left" + '\n')
                        # Right
                        elif row["x"] > x_vals_max:
                            tri_value_list.append(rgi((x_vals_max, row["y"], row["depth"])))
                            # f.write(str(row["x"]) + ", " + str(row["y"]) + ", " + str(
                            #     rgi((x_vals_max, row["y"], row["depth"]))) + ", Right" + '\n')
                    else:
                        # f.write(str(row["x"]) + ", " + str(row["y"]) + ", 100000" + ', Flag 5\n')
                        pass
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
            log_file("Completed block: " + str(input_bathymetry_split_i) + "\n", output_directory,
                     input_bathymetry_list, verbose_mode)

            del tri_value_list, tri, depth_array, depth_array_min, depth_array_max, input_bathymetry_split
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
    print "... Mosaicking using GDAL"
    arcpy.env.overwriteOutput = True
    if not arcpy.Exists(os.path.join(output_directory, "Outputs_GDAL")):
        os.makedirs(os.path.join(output_directory, "Outputs_GDAL"))
    counter = 0
    for raster_output in output_list_chunk:
        counter += 1
        arcpy.Mirror_management(in_raster=os.path.join(output_directory, "Outputs", raster_output),
                                out_raster=os.path.join(output_directory, "Outputs_GDAL", "sp" + str(counter) + ".tif"))
    gdal_location = r"c:\OSGeo4W64\bin\gdalwarp.exe --config GDAL_DATA C:\OSGeo4W64\share\gdal "
    gdal_files = os.path.join(output_directory, "Outputs_GDAL", "sp*.tif")
    gdal_output = os.path.join(output_directory, "tri_output.tif")
    with open(os.devnull, "w") as f:
        subprocess.call(gdal_location + gdal_files + " " + gdal_output, stdout=f)
    arcpy.AddMessage("... Mosaicking using GDAL complete")

def mpchunk(output_directory, input_bathymetry_cs, output_list_chunk):
    print "... Processing chunk: " + str(output_list_chunk[0])
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

    print "... Completed mosaicking chunk: " + str(output_list_chunk[0])

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

def log_file(text, output_directory, filename, verbose_mode):
    log_file_write = open(os.path.join(output_directory, "Outputs", str(filename) + ".log"), "a")
    log_file_write.write(str(text))
    log_file_write.close()
    if verbose_mode:
        print str(text)

if __name__ == '__main__':


    env_data = [
        [r"D:\GlobENV\1_Input_Environmental_Datasets\Regional\Mediterranean\dox\dox_max\Projected", "dox", "doxmax"],
        [r"D:\GlobENV\1_Input_Environmental_Datasets\Regional\Mediterranean\dox\dox_mean\Projected", "dox", "doxmean"],
        [r"D:\GlobENV\1_Input_Environmental_Datasets\Regional\Mediterranean\dox\dox_min\Projected", "dox", "doxmin"],
        [r"D:\GlobENV\1_Input_Environmental_Datasets\Regional\Mediterranean\dox\dox_std\Projected", "dox", "doxstd"]
        ]

    bathymetries = [[r"D:\GlobENV\2_Bathymetries_and_Outputs\Mediterranean\emodnet18_sr", "emod"],
                    [r"D:\GlobENV\2_Bathymetries_and_Outputs\Mediterranean\geb20med_sr", "gebco"]
                    ]

    output_directory_base = r"D:\GlobENV\2_Bathymetries_and_Outputs\Mediterranean"

    cpu_cores_used = "10"
    chunk_mode = 'gdal'  # gdal or arcpy - GDAL way faster, need to install GDAL binaries in standard location
    test_mode = False
    verbose_mode = True

    for b in bathymetries:
        for e in env_data:
                globenv(b[0], e[0], e[1], os.path.join(output_directory_base, b[1], e[2]), cpu_cores_used, chunk_mode, test_mode,
                        verbose_mode)


    #
    # env_data = [
    #     [r"D:\GlobENV\1_Input_Environmental_Datasets\world-ocean_atlas-2018\temperature\Projected", "t_an", "temp"],
    #     [r"D:\GlobENV\1_Input_Environmental_Datasets\world-ocean_atlas-2018\salinity\Projected", "s_an", "sal"],
    #     [r"D:\GlobENV\1_Input_Environmental_Datasets\world-ocean_atlas-2018\aoxu\Output\Projected", "a_an", "aoxu"],
    #     [r"D:\GlobENV\1_Input_Environmental_Datasets\world-ocean_atlas-2018\nitrate\Output\Projected", "n_an", "nit"],
    #     [r"D:\GlobENV\1_Input_Environmental_Datasets\world-ocean_atlas-2018\oxygen\Projected", "o_an", "diso2"],
    #     [r"D:\GlobENV\1_Input_Environmental_Datasets\world-ocean_atlas-2018\phosphate\Output\Projected", "p_an", "phos"],
    #     [r"D:\GlobENV\1_Input_Environmental_Datasets\world-ocean_atlas-2018\pos\Output\Projected", "o_an", "pos"],
    #     [r"D:\GlobENV\1_Input_Environmental_Datasets\world-ocean_atlas-2018\silicate\Output\Projected", "i_an", "sil"]
    #     ]
    #
    # bathymetries = [[r"D:\GlobENV\2_Bathymetries_and_Outputs\NEA\srtm_sr", "srtm"],
    #                 [r"D:\GlobENV\2_Bathymetries_and_Outputs\NEA\gebco_sr", "gebco"]
    #                 ]
    #
    # output_directory_base = r"D:\GlobENV\2_Bathymetries_and_Outputs\NEA"
    #
    # cpu_cores_used = "10"
    # chunk_mode = 'gdal'  # gdal or arcpy - GDAL way faster, need to install GDAL binaries in standard location
    # test_mode = False
    # verbose_mode = False
    #
    # for b in bathymetries:
    #     for e in env_data:
    #         if b[1] == 'srtm' and e[2] == 'temp':
    #             pass
    #         elif b[1] == 'srtm' and e[2] == 'sal':
    #             pass
    #         elif b[1] == 'srtm' and e[2] == 'aoxu':
    #             pass
    #         elif b[1] == 'srtm' and e[2] == 'nit':
    #             pass
    #         else:
    #             globenv(b[0], e[0], e[1], os.path.join(output_directory_base, b[1], e[2]), cpu_cores_used, chunk_mode, test_mode,
    #                     verbose_mode)