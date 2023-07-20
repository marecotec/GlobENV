##################################################################################################
# GlobENV - Data Preparer Bathymetry - GlobENV_Data_Preparer_Bathymetry.py
#
# Dr Andy Davies - University of Rhode Island, USA - github: marecotec
#
# This python code requires ArcGIS 10.4+.
#
##################################################################################################
# Copyright 2019 Dr Andy Davies
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

#  GlobENV Bathymetry Data Preparation

import arcpy
import os, csv
from arcpy.sa import *
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


# This function pre-reads data to assess whether or not this script
# can handle it, and prepare it for GlobENV

def data_desc_raster(input_bathymetry, status_print):
    try:
        if input_bathymetry is not None:

            input_bathymetry_prop = {}

            if status_print == True:
                arcpy.AddMessage("Attempting to read and describe raster file: " + str(input_bathymetry))

            input_bathymetry_desc = arcpy.Describe(input_bathymetry)

            if input_bathymetry_desc.dataType == 'RasterDataset':

                if status_print == True:
                    arcpy.AddMessage("-- Raster format is: " + str(input_bathymetry_desc.format))
                    arcpy.AddMessage("-- Raster CS type: " + str(input_bathymetry_desc.spatialReference.type))
                    arcpy.AddMessage("-- Raster CS name: " + str(input_bathymetry_desc.spatialReference.name))
                    arcpy.AddMessage(
                        "-- Raster xy units: " + str(input_bathymetry_desc.spatialReference.linearUnitName))

                input_bathymetry_prop["extent"] = input_bathymetry_desc.extent
                input_bathymetry_prop["format"] = input_bathymetry_desc.format
                input_bathymetry_prop["cs type"] = input_bathymetry_desc.spatialReference.type
                input_bathymetry_prop["cs name"] = input_bathymetry_desc.spatialReference.name
                input_bathymetry_prop["xy units"] = input_bathymetry_desc.spatialReference.linearUnitName

                input_bathymetry_prop["cell size x"] = arcpy.GetRasterProperties_management(input_bathymetry,
                                                                                            "CELLSIZEX").getOutput(0)
                input_bathymetry_prop["cell size y"] = arcpy.GetRasterProperties_management(input_bathymetry,
                                                                                            "CELLSIZEY").getOutput(0)
                input_bathymetry_prop["minimum"] = arcpy.GetRasterProperties_management(input_bathymetry,
                                                                                        "MINIMUM").getOutput(0)
                input_bathymetry_prop["maximum"] = arcpy.GetRasterProperties_management(input_bathymetry,
                                                                                        "MAXIMUM").getOutput(0)
                input_bathymetry_prop["nrows"] = arcpy.GetRasterProperties_management(input_bathymetry,
                                                                                      "ROWCOUNT").getOutput(0)
                input_bathymetry_prop["ncols"] = arcpy.GetRasterProperties_management(input_bathymetry,
                                                                                      "COLUMNCOUNT").getOutput(0)

                if input_bathymetry_prop.get("cell size x", "1") != input_bathymetry_prop.get("cell size y", "2"):
                    arcpy.AddError(
                        "Warning: Bathymetry x and y cell sizes should be identical (i.e. data must be in square cells).")
                    sys.exit(0)

                if status_print == True:
                    arcpy.AddMessage(
                        "-- Raster cell size: " + str(round(float(input_bathymetry_prop.get("cell size x", "!!!")), 2)))
                    arcpy.AddMessage("-- Raster rows: " + str(input_bathymetry_prop.get("nrows", "!!!")))
                    arcpy.AddMessage("-- Raster cols: " + str(input_bathymetry_prop.get("ncols", "!!!")))
                    arcpy.AddMessage("-- Estimated cell count: " + str(
                        int(input_bathymetry_prop.get("nrows", "!!!")) * int(input_bathymetry_prop.get("nrows",
                                                                                                       "none"))) + " .. does not take into account no_data cells")

                    if float(input_bathymetry_prop.get("minimum", "!!!")) < 0:
                        arcpy.AddMessage("-- Depth range: " + str(
                            round(float(input_bathymetry_prop.get("maximum", "!!!")), 2)) + " - " + str(
                            round(float(input_bathymetry_prop.get("minimum", "!!!")), 2)))
                    else:
                        arcpy.AddMessage("-- Depth range: " + str(
                            round(float(input_bathymetry_prop.get("minimum", "!!!")), 2)) + " - " + str(
                            round(float(input_bathymetry_prop.get("maximum", "!!!")),
                                  2)) + " .. you should use negative bathymetry if possible")

                validity = True

            else:
                arcpy.AddError("Bathymetry data is not in raster format")
                validity = False

        else:
            arcpy.AddMessage(".. file read unsuccessful, ")
            validity = False
    except:
        validity = False

    return input_bathymetry_prop, validity


def process_bathymetry(input_bathymetry, output_directory, test_mode):
    arcpy.AddMessage("GlobENV - Bathymetry data preparation script..")

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    if not os.path.exists(os.path.join(output_directory, "temp")):
        os.makedirs(os.path.join(output_directory, "temp"))

    arcpy.AddMessage(".. preparing data, if your raster is large, this step will take a short while...")

    if not arcpy.Exists(os.path.join(output_directory, "temp", "rst")):
        arcpy.CopyRaster_management(input_bathymetry, os.path.join(output_directory, "temp", "rst"), format="Esri Grid")

    input_bathymetry = os.path.join(output_directory, "temp", "rst")

    input_bathymetry_prop, validity = data_desc_raster(input_bathymetry, True)

    arcpy.AddMessage("Now checking raster file: " + str(input_bathymetry))

    if input_bathymetry_prop.get("cs type", "!!!") == "Projected":
        arcpy.AddMessage("-- Is projected: " + "Y")
    else:
        arcpy.AddMessage("-- Is projected: " + "N")
        sys.exit(0)

    if input_bathymetry_prop.get("xy units", "!!!") == "Meter":
        arcpy.AddMessage("-- Is in meters: " + "Y")
    else:
        arcpy.AddMessage("-- Is in meters: " + "N")
        sys.exit(0)

    if validity:
        arcpy.AddMessage("-- Raster is valid: " + "Y")
        input_bathymetry = os.path.join(output_directory, "temp", "rst")
    else:
        arcpy.AddMessage("-- Raster is valid: " + "N")
        sys.exit(0)

    # If the clauses above are met, we now move onto processing

    if int(input_bathymetry_prop.get("nrows", "!!!")) < 1999 & int(input_bathymetry_prop.get("ncols", "!!!")) < 1999:
        # In this case, only a single raster needs to be created as it falls below our depth block size threshold for
        # memory limitations in ArcPy 32bit
        arcpy.AddMessage(
            "Now processing raster file... You will see some warnings about BUILDVAT, these can be ignored.")

        if not arcpy.Exists(os.path.join(output_directory, "sp0.asc")):
            arcpy.RasterToASCII_conversion(input_bathymetry, os.path.join(output_directory, "sp0.asc"))
            try:
                input_bathymetry_test = data_desc_raster(os.path.join(output_directory, "sp0.asc"), False)
                arcpy.AddMessage("Now testing raster file... sp0.asc" + " --- pass")
            except:
                arcpy.AddMessage("Now testing raster file... sp0.asc" + " --- fail")

        arcpy.CopyRaster_management(os.path.join(output_directory, input_bathymetry),
                                    os.path.join(output_directory, "bath_config"))

        arcpy.CopyRaster_management(os.path.join(output_directory, input_bathymetry),
                                    os.path.join(output_directory, "sp0"))


    else:
        # In this case, we need to split our raster into blocks
        arcpy.AddMessage("Now processing raster file...")

        if not arcpy.Exists(os.path.join(output_directory, "anull_rst")):
            arcpy.gp.SetNull_sa(input_bathymetry, input_bathymetry, os.path.join(output_directory, "anull_rst"),
                                "VALUE >= 0.0")
            input_bathymetry = os.path.join(output_directory, "anull_rst")

        if not arcpy.Exists(os.path.join(output_directory, "null_rst")):
            null_raster = Con(IsNull(input_bathymetry), 3000000, input_bathymetry)
            null_raster.save(os.path.join(output_directory, "null_rst.tif"))

        # Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
        # The following inputs are layers or table views: "sp12", "sp12"

        input_bathymetry = os.path.join(output_directory, "null_rst.tif")

        rExt = input_bathymetry_prop.get("extent", "!!!")

        rWid = rExt.XMax - rExt.XMin
        rHgt = rExt.YMax - rExt.YMin

        # try to get step size set correctly for both x and y
        nrows = int(input_bathymetry_prop.get("nrows", "!!!"))
        ncols = int(input_bathymetry_prop.get("ncols", "!!!"))

        # try to get cell size set correctly for both x and y
        rowcs = float(input_bathymetry_prop.get("cell size y", "!!!"))
        colcs = float(input_bathymetry_prop.get("cell size x", "!!!"))

        def constrainsteps(n):
            total = 1
            while n >= 1999:
                n = (n // 2)
                total += 2
            total = total * 2
            return total, n

        nrow_steps, nrow_n = constrainsteps(nrows)
        ncol_steps, ncol_n = constrainsteps(ncols)

        # arcpy.AddMessage("Resulting bathymetric tiles will be: " + str(nrow_n) + " by " + str(ncol_n) + " pixels.")

        xStep = rWid / ncol_steps  # the value of each step
        yStep = rHgt / nrow_steps  # in the X and Y

        while xStep / colcs > 4000.0:
            xStep = xStep / 2
            ncol_steps = ncol_steps * 2

        while yStep / rowcs > 4000.0:
            yStep = yStep / 2
            nrow_steps = nrow_steps * 2

        count = 0

        arcpy.CreateFeatureclass_management(os.path.join(output_directory), "polygon.shp", geometry_type="Polygon")
        arcpy.AddField_management(os.path.join(output_directory, "polygon.shp"), "Name", "TEXT", "", "", "25", "",
                                  "NON_NULLABLE", "NON_REQUIRED", "")

        for xIndex in range(ncol_steps):  # python iterating for var in list:
            for yIndex in range(nrow_steps):
                # calculate the extent, if you want overlap you can do it here
                Xmin = rExt.XMin + (xStep * xIndex)
                Xmax = Xmin + xStep
                Ymin = rExt.YMin + (yStep * yIndex)
                Ymax = Ymin + yStep

                # make the name of the out file using % formatting.
                OutName = "sp" + str(count)
                # define the clipping box
                ClipBox = "%f %f %f %f" % (Xmin, Ymin, Xmax, Ymax)

                pol_array = arcpy.Array([arcpy.Point(Xmin, Ymin),
                                         arcpy.Point(Xmin, Ymax),
                                         arcpy.Point(Xmax, Ymax),
                                         arcpy.Point(Xmax, Ymin)])
                polygon = arcpy.Polygon(pol_array)

                # do the clip using the box and out name into the out folder
                if test_mode:
                    arcpy.AddMessage("Creating: " + str(OutName) + ", with bounds: " + str(ClipBox))

                if not test_mode:
                    arcpy.Clip_management(input_bathymetry, ClipBox, os.path.join(output_directory, OutName))

                cursor = arcpy.da.InsertCursor(os.path.join(output_directory, "polygon.shp"), ['SHAPE@', 'Name'])
                cursor.insertRow([polygon, OutName])
                del cursor

                count = count + 1

        if not test_mode:
            arcpy.Delete_management(os.path.join(output_directory, "null_rst.tif"))
            arcpy.Delete_management(os.path.join(output_directory, "anull_rst"))

        env.workspace = output_directory

        if not test_mode:
            input_bathymetry_list = arcpy.ListRasters("*", "ALL")

            arcpy.AddMessage("There are " + str(len(input_bathymetry_list)) + " depth blocks created.")
            check_list = []
            pass_list = []
            null_list = []

            bath_config = False

            for i in input_bathymetry_list:

                min_test = float(arcpy.GetRasterProperties_management(i, "MINIMUM").getOutput(0)) + float(
                    arcpy.GetRasterProperties_management(i, "MAXIMUM").getOutput(0))

                if min_test > float(5999999):
                    arcpy.AddMessage("Now testing raster file..." + i + "-- all null")
                    arcpy.Delete_management(os.path.join(output_directory, i))
                    null_list.append(i)

                else:
                    arcpy.gp.SetNull_sa(i, i, os.path.join(output_directory, i + "_1"),
                                        "VALUE = 3000000")

                    arcpy.Delete_management(os.path.join(output_directory, i))
                    arcpy.Rename_management(os.path.join(output_directory, i + "_1"), os.path.join(output_directory, i))
                    arcpy.Delete_management(os.path.join(output_directory, i + "_1"))

                    # 1 test depth block to make sure data is in it
                    try:
                        input_bathymetry_test, validity = data_desc_raster(i, False)

                        if validity:
                            arcpy.AddMessage("Now testing raster file..." + i + "-- pass")
                            pass_list.append(i)

                            if not bath_config:
                                arcpy.CopyRaster_management(os.path.join(output_directory, i),
                                                            os.path.join(output_directory, "bath_config"))
                                bath_config = True
                    except:
                        arcpy.AddMessage("Now testing raster file..." + i + "-- fail")
                        check_list.append(i)

            if len(check_list) > 0:
                arcpy.AddMessage(str(len(
                    check_list)) + " rasters failed the basic check, this likely means that they are invalid in terms of bathymetric data, " +
                                 "you should check these rasters for errors. They have been moved to the *output_directory*/FailedCheck directory.")

                arcpy.AddMessage("For reference, the following rasters failed checks: " + str(', '.join(check_list)))

                if not os.path.exists(os.path.join(output_directory, "FailedCheck")):
                    os.makedirs(os.path.join(output_directory, "FailedCheck"))

                for i in check_list:
                    arcpy.Copy_management(os.path.join(output_directory, i),
                                          os.path.join(output_directory, "FailedCheck", i), )
                    arcpy.Delete_management(os.path.join(output_directory, i))

                if len(null_list) > 0:
                    arcpy.AddMessage(
                        "For reference, the following rasters were completely null: " + str(', '.join(null_list)))

                arcpy.AddMessage(
                    "For reference, " + str(len(pass_list)) + " rasters passed tests and are considered valid")
                arcpy.AddMessage("Process complete with errors..")

            else:
                if len(null_list) > 0:
                    arcpy.AddMessage(
                        "For reference, the following rasters were completely null: " + str(', '.join(null_list)))

                arcpy.AddMessage(
                    "For reference, " + str(len(pass_list)) + " rasters passed tests and are considered valid")
                arcpy.AddMessage("Process complete without errors..")
        else:
            arcpy.AddMessage("Test mode complete... Check your shapefile")


if __name__ == '__main__':

    test_mode = False

    input_bathymetry = r"G:\GEBCO_2022\gebco_2022_ocean_Merc.tif"
    output_directory = r"G:\GlobENV_SRTM15\Processing\2_Bathymetries_and_Outputs\gebco2022\gebco_sr"

    process_bathymetry(input_bathymetry, output_directory, test_mode)
