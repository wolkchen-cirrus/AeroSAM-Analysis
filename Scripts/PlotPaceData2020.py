"""
PlotPaceData2019 - This script will use the AirborneParticleAnalysis library to perform the level 2 analysis on the PaCE
2019 data in particular. pacePlots is the module containing the plotting functions, level1to2 is the module containing
the analysis functions and objects for deriving level 2 data, and common contains all the misc functions and objects
shared between multiple modules/scripts. Most of the body of this script is an iterative searching algorithm for two
different data sets. Different algorithms were considered, but brevity was also.
"""

from AirborneParticleAnalysis import PacePlots
from AirborneParticleAnalysis import level1to2
from AirborneParticleAnalysis import common
import numpy as np
from os import listdir
from matplotlib import pyplot as plt


if __name__ == "__main__":
    print "============ Executing the Plotting Script for the PaCE 2020 Campaign Data ============\n"

    # Starting the import of the level 2 data, the script "convert_all_0to1.py" must be run first to obtain this level 1
    # data in the first instance.
    data_dir = common.read_setting("base_data_dir")
    t1 = listdir(data_dir)
    data_files_talon = []
    for date in t1:
        date_dir = data_dir + "\\" + date + "\\" + "level_2"
        t2 = listdir(date_dir)
        for fnm in t2:
            if "FMITalon_" in fnm:
                data_files_talon.append(date_dir + "\\" + fnm)
    data_files_static = []
    for date in t1:
        date_dir = data_dir + "\\" + date + "\\" + "level_1"
        t2 = listdir(date_dir)
        for fnm in t2:
            if "StaticUCASS_" in fnm:
                data_files_static.append(date_dir + "\\" + fnm)

    data_dict_talon = {}
    for file_name in data_files_talon:
        print "INFO: Importing file - %s" % file_name
        key_name = file_name.split("\\")[-1].split(".")[0]  # Dict key name is file name
        data_dict_talon[key_name] = level1to2.import_level1(file_name)  # Importing, uses "cPickle"
    data_dict_static = {}
    for file_name in data_files_static:
        print "INFO: Importing file - %s" % file_name
        key_name = file_name.split("\\")[-1].split(".")[0]  # Dict key name is file name
        data_dict_static[key_name] = level1to2.import_level1(file_name)  # Importing, uses "cPickle"

    station_altitude_asl_m = float(common.read_setting("station_altitude_asl_mm")) / 1000.0

    talon_bbs = None
    static_bbs = None
    sorted_t_dn_arr = []
    sorted_s_dn_arr = []
    sorted_dt_arr = []
    for data_file_talon in data_dict_talon:

        # Matching data files together
        date_time_talon = str(data_file_talon.split("_")[-3]) + str(data_file_talon.split("_")[-2])
        diff_arr = []
        for data_file_static in data_dict_static:
            date_time_static = str(data_file_static.split("_")[-3]) + str(data_file_static.split("_")[-2])
            diff_arr.append(int(date_time_talon) - int(date_time_static))
        min_index = diff_arr.index(min(i for i in diff_arr if i > 0))
        matched_static_file = data_dict_static[data_dict_static.keys()[min_index]]

        sequence = ["Up", "Down"]
        prof_num = data_dict_talon[data_file_talon].up_profile_mask.shape[1]

        for pn in range(prof_num):
            for j in sequence:

                date_time_string = str(data_file_talon.split("_")[-3])[:4] + "/" + \
                                   str(data_file_talon.split("_")[-3])[4:6] + "/" + \
                                   str(data_file_talon.split("_")[-3])[6:] + " " + \
                                   str(data_file_talon.split("_")[-2])[:2] + ":" + \
                                   str(data_file_talon.split("_")[-2])[2:4]
                sorted_dt_arr.append(date_time_string)

                # Getting matching altitude index for Talon data
                row_number = level1to2.fetch_row(altitude=station_altitude_asl_m,
                                                 level1_data=data_dict_talon[data_file_talon],
                                                 profile=j, prof_num=pn)[0]
                rows_for_mean_talon = level1to2.fetch_row_tolerance(altitude=station_altitude_asl_m,
                                                                    level1_data=data_dict_talon[data_file_talon],
                                                                    profile=j, prof_num=pn)

                # Getting time the the SUA was at station altitude
                sync_time = level1to2.get_time_from_alt(data_dict_talon[data_file_talon], row_number)

                # Getting matching time index for static UCASS data
                rows_for_mean_static = level1to2.fetch_row_tolerance(time=sync_time, level1_data=matched_static_file)

                # Calculating properties from located times and altitudes (Talon)
                mean_dn_dlogdp_talon = level1to2.mean_dn_dlogdp(data_dict_talon[data_file_talon], rows_for_mean_talon)
                row_index = [np.where(np.asarray(data_dict_talon[data_file_talon].alt) == i)
                             [0][0] for i in [k[0] for k in rows_for_mean_talon]]
                n_conc_buf_talon = np.true_divide(np.mean(data_dict_talon[data_file_talon]
                                                          .number_concentration[row_index]), 1e6)

                # Calculating properties from located times and altitudes (Static)
                mean_dn_dlogdp_static = level1to2.mean_dn_dlogdp(matched_static_file, rows_for_mean_static)
                row_index = [np.where(np.asarray(matched_static_file.time) == i)
                             [0][0] for i in [k[0] for k in rows_for_mean_static]]
                n_conc_buf_static = np.true_divide(np.mean(matched_static_file.number_concentration[row_index]), 1e6)

                talon_bbs = data_dict_talon[data_file_talon].bin_bounds_dp_um
                static_bbs = matched_static_file.bin_bounds_dp_um

                sorted_t_dn_arr.append(mean_dn_dlogdp_talon)
                sorted_s_dn_arr.append(mean_dn_dlogdp_static)

    PacePlots.plot_pace_dn_dlogdp_2020(sorted_t_dn_arr, sorted_s_dn_arr, talon_bbs, static_bbs, sorted_dt_arr)
    plt.show()
    pass
