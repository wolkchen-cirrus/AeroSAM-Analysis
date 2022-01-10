"""
PlotPaceData2020 - This script will use the AirborneParticleAnalysis library to perform the level 2 analysis on the PaCE
2020 data in particular. pacePlots is the module containing the plotting functions, level1to2 is the module containing
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
import datetime
import os


def _get_wind_direction(level1_object):
    filename_date = level1_object.path.split("\\")[-1].split("_")[-3]
    d_dir = common.read_setting("Sammal_Metd_Path")
    d_files = os.listdir(d_dir)
    date_diff = []
    for d_file in d_files:
        d_date = str(d_file.split("_")[-1].split(".")[0])
        date_diff.append(abs(float(d_date) - float(filename_date)))
    val = min(date_diff)
    matched_date_index = date_diff.index(val)
    d_path = d_dir + "\\" + d_files[matched_date_index]

    met_time_col = common.fetch_column(d_path, 0, remove_r1=True)
    met_wind_col = common.fetch_column(d_path, 7, remove_r1=True)
    met_epoch_col = [(datetime.datetime.strptime(d, '%Y-%m-%d %H:%M:%S')
                      - datetime.datetime(1970, 1, 1) - datetime.timedelta(hours=3)).total_seconds()
                     for d in met_time_col]
    flight_epoch = level1_object.time[0]
    i1, _ = common.sync_data_point(flight_epoch, met_epoch_col)
    w_dir = float(met_wind_col[i1])
    return w_dir


if __name__ == "__main__":
    print "============ Executing the Plotting Script for the PaCE 2020 Campaign Data ============\n"

    # Starting the import of the level 2 data, the script "convert_all_0to1.py" must be run first to obtain this level 1
    # data in the first instance. "convert_FMI_1to2.py" must also be run to obtain level 2 FMI-Talon data.
    data_dir = common.read_setting("base_data_dir")
    station_altitude_asl_m = float(common.read_setting("station_altitude_asl_mm")) / 1000.0

    # Change these values to true to generate different plots.
    plot_d_eff = False
    plot_dn_dlog_dp = False
    plot_n_conc = True
    plot_rms_pry = True

    t1 = listdir(data_dir)
    d_eff_arr = []
    diff_n_conc_arr = []
    aoa_arr = []
    asp_arr = []
    rms_pry = []
    dt_arr = []
    for date in t1:
        data_files_talon = []
        data_talon = []
        date_dir = data_dir + "\\" + date + "\\" + "level_2"
        t2 = listdir(date_dir)
        for fnm in t2:
            if "FMITalon_" in fnm:
                data_talon.append(level1to2.import_level1(date_dir + "\\" + fnm))
        if not data_talon:
            continue
        data_files_static = []
        data_static = []
        date_dir = data_dir + "\\" + date + "\\" + "level_1"
        t2 = listdir(date_dir)
        for fnm in t2:
            if "StaticUCASS_" in fnm:
                data_static.append(level1to2.import_level1(date_dir + "\\" + fnm))
        if not data_static:
            continue

        talon_bbs = None
        static_bbs = None

        sorted_t_dn_arr = []
        sorted_s_dn_arr = []
        sorted_t_eff_arr = []
        sorted_s_eff_arr = []
        sorted_dt_arr = []
        for data_file_talon in data_talon:

            # Matching data files together
            date_time_talon = str(data_file_talon.path.split("_")[-3]) + str(data_file_talon.path.split("_")[-2])
            diff_arr = []
            for data_file_static in data_static:
                date_time_static = str(data_file_static.path.split("_")[-3]) + str(data_file_static.path.split("_")[-2])
                diff_arr.append(int(date_time_talon) - int(date_time_static))
            min_index = diff_arr.index(min(i for i in diff_arr if i > 0))
            matched_static_file = data_static[min_index]

            wd = _get_wind_direction(data_file_talon)
            _, _, aoa_ar = level1to2.check_valid_fixedwing(data_file_talon, wd,
                                                           airspeed_type="adjusted", aoa_lim_deg=10, airspeed_lim_ms=20)

            sequence = ["Up", "Down"]
            prof_num = data_file_talon.up_profile_mask.shape[1]

            for pn in range(prof_num):

                for j in sequence:

                    if j == "Up":
                        prof_mask = data_file_talon.up_profile_mask
                    elif j == "Down":
                        prof_mask = data_file_talon.down_profile_mask
                    else:
                        raise ValueError

                    t_centre = data_file_talon.bin_centres_dp_um
                    s_centre = matched_static_file.bin_centres_dp_um

                    date_time_string = str(data_file_talon.path.split("_")[-3])[:4] + "/" + \
                        str(data_file_talon.path.split("_")[-3])[4:6] + "/" + \
                        str(data_file_talon.path.split("_")[-3])[6:] + " " + \
                        str(data_file_talon.path.split("_")[-2])[:2] + ":" + \
                        str(data_file_talon.path.split("_")[-2])[2:4]
                    sorted_dt_arr.append(date_time_string)
                    dt_arr.append(date_time_string)

                    # Getting matching altitude index for Talon data
                    row_number = level1to2.fetch_row(altitude=station_altitude_asl_m,
                                                     level1_data=data_file_talon,
                                                     profile=j, prof_num=pn)[0]
                    rows_for_mean_talon = level1to2.fetch_row_tolerance(altitude=station_altitude_asl_m,
                                                                        level1_data=data_file_talon,
                                                                        profile=j, prof_num=pn)

                    # Getting time the the SUA was at station altitude
                    sync_time = level1to2.get_time_from_alt(data_file_talon, row_number)

                    # Getting matching time index for static UCASS data
                    rows_for_mean_static = level1to2.fetch_row_tolerance(time=sync_time,
                                                                         level1_data=matched_static_file)

                    # Calculating properties from located times and altitudes (Talon)
                    mean_dn_dlogdp_talon = level1to2.mean_dn_dlogdp(data_file_talon, rows_for_mean_talon)
                    row_index = [np.where(np.asarray(data_file_talon.alt) == i)
                                 [0][0] for i in [k[0] for k in rows_for_mean_talon]]
                    n_conc_buf_talon = np.true_divide(np.mean(data_file_talon
                                                              .number_concentration[row_index]), 1e6)

                    # Calculating properties from located times and altitudes (Static)
                    mean_dn_dlogdp_static = level1to2.mean_dn_dlogdp(matched_static_file, rows_for_mean_static)
                    row_index = [np.where(np.asarray(matched_static_file.time) == i)
                                 [0][0] for i in [k[0] for k in rows_for_mean_static]]
                    n_conc_buf_static = np.true_divide(np.mean(matched_static_file.number_concentration[row_index]),
                                                       1e6)

                    talon_bbs = data_file_talon.bin_bounds_dp_um
                    static_bbs = matched_static_file.bin_bounds_dp_um

                    talon_bcs = data_file_talon.bin_centres_dp_um
                    static_bcs = matched_static_file.bin_centres_dp_um

                    mean_counts_talon = np.mean(
                        data_file_talon.raw_counts[[j for (i, j) in rows_for_mean_talon], :], axis=0)
                    mean_counts_static = np.mean(
                        matched_static_file.raw_counts[[j for (i, j) in rows_for_mean_static], :], axis=0)
                    mean_sv_talon = np.mean(data_file_talon.sample_volume_m3
                                            [[j for (i, j) in rows_for_mean_talon], 0], axis=0)
                    mean_sv_static = np.mean(matched_static_file.sample_volume_m3
                                             [[j for (i, j) in rows_for_mean_static], 0], axis=0)
                    mean_n_conc_talon = np.mean(
                        data_file_talon.number_concentration[[j for (i, j) in rows_for_mean_talon], :], axis=0)
                    mean_n_conc_static = np.mean(
                        matched_static_file.number_concentration[[j for (i, j) in rows_for_mean_static], :], axis=0)
                    mean_aoa_talon = np.mean(aoa_ar[[j for (i, j) in rows_for_mean_talon], :], axis=0)
                    mean_asp_talon = np.mean(
                        data_file_talon.vz_cms[[j for (i, j) in rows_for_mean_talon], :], axis=0) / 100.0

                    sorted_t_eff_arr.append(level1to2.effective_diameter(mean_counts_talon, mean_sv_talon, t_centre))
                    sorted_s_eff_arr.append(level1to2.effective_diameter(mean_counts_static, mean_sv_static, s_centre))

                    d_eff_arr.append(level1to2.effective_diameter(mean_counts_talon, mean_sv_talon, t_centre) -
                                     level1to2.effective_diameter(mean_counts_static, mean_sv_static, s_centre))
                    diff_n_conc_arr.append(mean_n_conc_talon - mean_n_conc_static)
                    aoa_arr.append(mean_aoa_talon)
                    asp_arr.append(mean_asp_talon)

                    sorted_t_dn_arr.append(mean_dn_dlogdp_talon)
                    sorted_s_dn_arr.append(mean_dn_dlogdp_static)

                    rms_p = level1to2.noise_measure(data_file_talon.pitch[[j for (i, j) in rows_for_mean_talon], :], 20)
                    rms_r = level1to2.noise_measure(data_file_talon.roll[[j for (i, j) in rows_for_mean_talon], :], 20)
                    rms_y = level1to2.noise_measure(data_file_talon.yaw[[j for (i, j) in rows_for_mean_talon], :], 20)

                    rms_pry.append(rms_p * rms_r * rms_y)

        if plot_dn_dlog_dp is True:
            PacePlots.plot_pace_dn_dlogdp_2020(sorted_t_dn_arr, sorted_s_dn_arr, talon_bbs, static_bbs, sorted_dt_arr,
                                               sorted_s_eff_arr, sorted_t_eff_arr, talon_bcs, static_bcs)

    if plot_d_eff is True:
        PacePlots.eff_dia_plot_2020(d_eff_arr, [i.split(" ")[1] for i in dt_arr])

    if plot_n_conc is True:
        PacePlots.conc_asp_2020([abs(i[0]/1000000.0) for i in diff_n_conc_arr[:8]], [i[0] for i in asp_arr[:8]],
                                [])

    if plot_rms_pry is True:
        PacePlots.conc_pry_2020([abs(i[0]/1000000.0) for i in diff_n_conc_arr][:], rms_pry[:],
                                [])

    plt.show()
    pass
