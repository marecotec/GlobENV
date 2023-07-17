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

def load_environmental_layers():

def load_bathymetry_layer():

def globenv_check():


if __name__ == '__main__':

    test_mode = False

    input_env_directory = r""
    input_bathymetry = r"M:\GlobEnv\2_Bathymetries_and_Outputs\Global\Gebco2020\gebco_2020_proj.tif"

    process_bathymetry(input_bathymetry, output_directory, test_mode)


