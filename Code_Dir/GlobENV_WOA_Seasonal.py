from GlobENV_64_bit import globenv
import os
import arcpy

#WOA Seasonal
env_data = [[r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\aoxu\autumn", "a_an", "aoxu_autumn"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\aoxu\spring", "a_an", "aoxu_spring"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\aoxu\summer", "a_an", "aoxu_summer"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\aoxu\winter", "a_an", "aoxu_winter"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\diso2\autumn", "o_an", "diso2_autumn"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\diso2\spring", "o_an", "diso2_spring"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\diso2\summer", "o_an", "diso2_summer"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\diso2\winter", "o_an", "diso2_winter"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\nit\autumn", "n_an", "nit_autumn"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\nit\spring", "n_an", "nit_spring"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\nit\summer", "n_an", "nit_summer"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\nit\winter", "n_an", "nit_winter"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\phos\autumn", "p_an", "phos_autumn"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\phos\spring", "p_an", "phos_spring"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\phos\summer", "p_an", "phos_summer"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\phos\winter", "p_an", "phos_winter"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\pos\autumn", "o_an", "pos_autumn"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\pos\spring", "o_an", "pos_spring"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\pos\summer", "o_an", "pos_summer"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\pos\winter", "o_an", "pos_winter"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\sal\autumn", "s_an", "sal_autumn"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\sal\spring", "s_an", "sal_spring"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\sal\summer", "s_an", "sal_summer"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\sal\winter", "s_an", "sal_winter"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\sil\autumn", "i_an", "sil_autumn"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\sil\spring", "i_an", "sil_spring"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\sil\summer", "i_an", "sil_summer"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\sil\winter", "i_an", "sil_winter"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\temp\autumn", "t_an", "temp_autumn"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\temp\spring", "t_an", "temp_spring"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\temp\summer", "t_an", "temp_summer"],
            [r"E:\1_Input_Environmental_Datasets_Mobile\world-ocean-atlas-2018-seasonal\temp\winter", "t_an", "temp_winter"]
            ]

bathymetries = [[r"E:\2_Bathymetries_and_Outputs\Other\North_Atlantic_SponGES_Matter\SR_Bathymetries\gebco_sr", "gebco"],
                [r"E:\2_Bathymetries_and_Outputs\Other\North_Atlantic_SponGES_Matter\SR_Bathymetries\srtm_sr", "srtm"]
                ]

output_directory_base = r"E:\2_Bathymetries_and_Outputs\Other\North_Atlantic_SponGES_Matter\WOASeasonal"

cpu_cores_used = "23"
chunk_mode = 'gdal'  # gdal or arcpy - GDAL way faster, need to install GDAL binaries in standard location
test_mode = False
verbose_mode = False

for b in bathymetries:
    for e in env_data:
            globenv(b[0], e[0], e[1], os.path.join(output_directory_base, b[1], e[2]), cpu_cores_used, chunk_mode, test_mode,
                    verbose_mode)
            if chunk_mode == 'gdal':
                arcpy.Copy_management(os.path.join(output_directory_base, b[1], e[2], "tri_output.tif"),
                                      os.path.join(output_directory_base, b[1], e[2] + ".tif"))