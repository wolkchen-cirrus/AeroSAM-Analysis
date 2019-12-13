import importer
import level0to1
import common
import os


if __name__ == "__main__":
    data_dir = common.read_setting("base_data_dir")
    data_dates = os.listdir(data_dir)
    if data_dir[-1] != "\\":
        data_dir += "\\"
    for date in data_dates:
        level0_path = data_dir + date + "\\" + "level_0" + "\\"
        level1_path = data_dir + date + "\\" + "level_1" + "\\"
        level0_files = os.listdir(level0_path)
        level1_files = os.listdir(level1_path)
        level1_checks = [i.split("_")[-2] for i in level1_files]
        convert = False
        convert_list = []
        for file_0 in level0_files:
            check_val = file_0.split("_")[-2]
            if check_val in level1_checks:
                convert_list.append(file_0)
            else:
                convert = True
        if convert is True:
            for file_0 in level0_files:
                data0_file_path = level0_path + file_0
                if file_0 not in convert_list:
                    if "CAS" in file_0:
                        level0_object = importer.StaticCASData(level0_path=data0_file_path)
                        level0to1.bin_centre_dp_um(level0_object)
                        level0to1.sample_volume(level0_object)
                        level0to1.dn_dlogdp(level0_object)
                        level0to1.export_level1(level0_object)
                    elif "AeroSAM-log" in file_0:
                        level0_object = importer.SUAData(level0_path=data0_file_path)
                        if level0_object.trash is True:
                            print "INFO: Trash data with filename %s, not converting to level 1" % file_0
                            continue
                        level0to1.split_by_pressure(level0_object)
                        level0to1.assign_ucass_lut(level0_object)
                        level0to1.bin_centre_dp_um(level0_object)
                        level0to1.sample_volume(level0_object)
                        level0to1.mass_concentration_kgm3(level0_object)
                        level0to1.num_concentration_m3(level0_object)
                        level0to1.dn_dlogdp(level0_object)
                        level0to1.export_level1(level0_object)
                    else:
                        raise RuntimeError("ERROR: Unrecognised data file")
                else:
                    print "INFO: Filename %s already converted" % file_0
                    continue
