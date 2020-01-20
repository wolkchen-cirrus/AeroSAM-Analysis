from AirborneParticleAnalysis import PacePlots
from AirborneParticleAnalysis import level1to2
from os import listdir


def _hhmmss_to_sec(hhmmss):
    return int(hhmmss[0:2])*3600+int(hhmmss[2:4])*60+int(hhmmss[4:6])


if __name__ == "__main__":
    print "============ Executing the Plotting Script for the PaCE 2019 Campaign Data ============\n"

    data_dir = "C:\\Users\\jg17acv\\University of Hertfordshire\\AeroSAM - Documents\\Data\\2019\\19-09-28\\level_1"
    station_altitude_asl_mm = 560000
    data_files = listdir(data_dir)
    data_dict = {}
    for file_name in data_files:
        print "INFO: Importing file - %s" % file_name
        key_name = file_name.split(".")[0]
        file_path = data_dir + "\\" + file_name
        data_dict[key_name] = level1to2.import_level1(file_path)

    fssp_bins = None
    cas_bins = None
    sam_bins = None
    dn_dlogdp_sam = {}
    dn_dlogdp_cas = {}
    dn_dlogdp_fssp = {}
    times = []
    print "INFO: Analysis pipeline first pass"
    for key in data_dict:
        if ("CAS" in key) or ("FSSP" in key):
            print "INFO: Skipping CAS data and FSSP data on first round"
            continue
        else:
            print "INFO: Processing SAM data"
            dn_buf = data_dict[key].dn_dlogdp
            row_number = level1to2.fetch_row(altitude=station_altitude_asl_mm, level1_data=data_dict[key])[0]
            times.append(_hhmmss_to_sec(data_dict[key].datetime[-6:]))
            dp_key = key + "_" + str(_hhmmss_to_sec(data_dict[key].datetime[-6:]))
            dn_dlogdp_sam[dp_key] = dn_buf[row_number]
            if sam_bins is None:
                sam_bins = data_dict[key].bin_centres_dp_um

    print "INFO: Analysis pipeline Second pass"
    for time in times:
        for key in data_dict:
            if ("CAS" in key) or ("FSSP" in key):
                print "INFO: Processing CAS and FSSP data"
                dn_buf = data_dict[key].dn_dlogdp
                row_number = level1to2.fetch_row(time=time, level1_data=data_dict[key])[0]
                dn_key = key + "_" + str(time)
                if "CAS" in key:
                    cas_bins = data_dict[key].bin_centres_dp_um
                    dn_dlogdp_cas[dn_key] = dn_buf[row_number]
                elif "FSSP" in key:
                    fssp_bins = data_dict[key].bin_centres_dp_um
                    dn_dlogdp_fssp[dn_key] = dn_buf[row_number]
            else:
                print "INFO: Skipping SAM data on second round"
                continue

    dn_dlogdp_comp = {}
    for sam_key in dn_dlogdp_sam:
        time = sam_key.split("_")[-1]
        buf = {sam_key: dn_dlogdp_sam[sam_key]}
        for cas_key in dn_dlogdp_cas:
            if time == cas_key.split("_")[-1]:
                buf[cas_key] = dn_dlogdp_cas[cas_key]
                break
        for fssp_key in dn_dlogdp_fssp:
            if time == fssp_key.split("_")[-1]:
                buf[fssp_key] = dn_dlogdp_fssp[fssp_key]
                break
        dn_dlogdp_comp[time] = buf

    for time in dn_dlogdp_comp:
        PacePlots.plot_pace_dn_dlogdp(dn_dlogdp_comp[time], sam_bins=sam_bins, cas_bins=cas_bins, fssp_bins=fssp_bins)
    pass
