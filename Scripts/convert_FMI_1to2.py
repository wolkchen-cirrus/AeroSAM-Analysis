from AirborneParticleAnalysis import level1to2, level0to1, common
import numpy as np
import datetime
import os


if __name__ == "__main__":
    data_dir = common.read_setting("base_data_dir")
    data_dates = os.listdir(data_dir)
    if data_dir[-1] != "\\":
        data_dir += "\\"
    for date in data_dates:
        level1_path = data_dir + date + "\\" + "level_1" + "\\"
        level2_path = data_dir + date + "\\" + "level_2" + "\\"
        level1_files = os.listdir(level1_path)
        level2_files = os.listdir(level2_path)
        level2_checks = ["_".join([i.split("_")[0], i.split("_")[-2]]) for i in level2_files]
        convert = False
        convert_list = []
        for file_0 in level1_files:
            check_val = "_".join([file_0.split("_")[0], file_0.split("_")[-2]])
            if check_val in level2_checks:
                convert_list.append(file_0)
            else:
                convert = True
        if convert is True:
            for file_1 in level1_files:
                data1_file_path = level1_path + file_1
                if file_1 not in convert_list:
                    if "FMITalon_" in file_1:
                        level1_object = level1to2.import_level1(data1_file_path)

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
                        wd = float(met_wind_col[i1])

                        level1to2.adjust_all_airspeed_mtof(level1_object)
                        level1to2.check_valid_fixedwing(level1_object, wd, airspeed_type="adjusted", aoa_lim_deg=15,
                                                        airspeed_lim_ms=20)

                        level0to1.sample_volume(level1_object, airspeed_type="adjusted")
                        level0to1.mass_concentration_kgm3(level1_object)
                        level0to1.num_concentration_m3(level1_object)
                        level0to1.dn_dlogdp(level1_object)
                        level1to2.export_level2(level1_object)

                    elif "StaticUCASS" in file_1:
                        level1_object = level1to2.import_level1(data1_file_path)

                        # ToDo: Realised that this has to be implemented in importer as part of level 1 analysis since
                        #  the column object doesn't work this way.
                        # Get mask of all rows to remove
                        period = level1_object.opc_aux[:, 0]
                        valid_mask = np.zeros(period.shape)
                        for p, i in zip(period, range(period.shape[0])):
                            if (int(p[0]) == 0) or (int(p[0]) == 255):
                                valid_mask[i] = 1
                                try:
                                    valid_mask[i+1] = 1
                                except IndexError:
                                    pass
                            else:
                                pass

                        # Count through arrays backwards to remove indexing errors
                        for v, i in zip(np.flip(valid_mask, 0), range(valid_mask.shape[0]-1, -1, -1)):
                            if v == 1:
                                pass
                            else:
                                pass

                    else:
                        print ("INFO: Skipping data file of unknown type")
                else:
                    print "INFO: Filename %s already converted" % file_1
                    continue
