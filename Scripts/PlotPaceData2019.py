from AirborneParticleAnalysis import PacePlots
from AirborneParticleAnalysis import level1to2
from AirborneParticleAnalysis import common
from os import listdir


plot_mean_dn_dlogdp = True
plot_exact_dn_dlogdp = False


def _hhmmss_to_sec(hhmmss):
    return int(hhmmss[0:2])*3600+int(hhmmss[2:4])*60+int(hhmmss[4:6])


if __name__ == "__main__":
    print "============ Executing the Plotting Script for the PaCE 2019 Campaign Data ============\n"

    data_dir = "C:\\Users\\jg17acv\\University of Hertfordshire\\AeroSAM - Documents\\Data\\2019\\19-09-28\\level_1"
    station_altitude_asl_mm = float(common.read_setting("station_altitude_asl_mm"))
    data_files = listdir(data_dir)
    data_dict = {}
    for file_name in data_files:
        print "INFO: Importing file - %s" % file_name
        key_name = file_name.split(".")[0]
        file_path = data_dir + "\\" + file_name
        data_dict[key_name] = level1to2.import_level1(file_path)

    fssp_bins = None
    cas_bins = None
    sam_bins_droplet = None
    sam_bins_aerosol = None
    dn_dlogdp_sam = {}
    dn_dlogdp_cas = {}
    dn_dlogdp_fssp = {}
    dn_dlogdp_sam_mean = {}
    dn_dlogdp_cas_mean = {}
    dn_dlogdp_fssp_mean = {}
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
            rows_for_mean = level1to2.fetch_row_tolerance(altitude=station_altitude_asl_mm, level1_data=data_dict[key])
            mean_dn_dlogdp = level1to2.mean_dn_dlogdp(data_dict[key], rows_for_mean[0])

            current_time = level1to2.get_time_from_alt(data_dict[key], row_number) + 3600   # 1hour time difference
            times.append(current_time)

            tags = data_dict[key].tags
            if "Droplet" in tags:
                sam_bins_droplet = data_dict[key].bin_centres_dp_um
                dn_key = key + "_" + "Droplet"
            elif "Aerosol" in tags:
                sam_bins_aerosol = data_dict[key].bin_centres_dp_um
                dn_key = key + "_" + "Aerosol"
            else:
                raise AttributeError("ERROR: Gain mode not specified")

            dn_key = dn_key + "_" + str(current_time)
            dn_dlogdp_sam[dn_key] = dn_buf[row_number]
            dn_dlogdp_sam_mean[dn_key] = mean_dn_dlogdp

    print "INFO: Analysis pipeline Second pass"
    for time in times:
        for key in data_dict:
            if ("CAS" in key) or ("FSSP" in key):
                print "INFO: Processing CAS and FSSP data"
                dn_buf = data_dict[key].dn_dlogdp
                row_number = level1to2.fetch_row(time=time, level1_data=data_dict[key])[0]
                rows_for_mean = level1to2.fetch_row_tolerance(time=time, level1_data=data_dict[key])
                mean_dn_dlogdp = level1to2.mean_dn_dlogdp(data_dict[key], rows_for_mean[0])
                dn_key = key + "_" + str(time)
                if "CAS" in key:
                    cas_bins = data_dict[key].bin_centres_dp_um
                    dn_dlogdp_cas[dn_key] = dn_buf[row_number]
                    dn_dlogdp_cas_mean[dn_key] = mean_dn_dlogdp
                elif "FSSP" in key:
                    fssp_bins = data_dict[key].bin_centres_dp_um
                    dn_dlogdp_fssp[dn_key] = dn_buf[row_number]
                    dn_dlogdp_fssp_mean[dn_key] = mean_dn_dlogdp
            else:
                print "INFO: Skipping SAM data on second round"
                continue

    print "INFO: Analysis pipeline third pass"
    dn_dlogdp_comp = {}
    dn_dlogdp_comp_mean = {}
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
    for sam_key in dn_dlogdp_sam_mean:
        time = sam_key.split("_")[-1]
        buf = {sam_key: dn_dlogdp_sam_mean[sam_key]}
        for cas_key in dn_dlogdp_cas_mean:
            if time == cas_key.split("_")[-1]:
                buf[cas_key] = dn_dlogdp_cas_mean[cas_key]
                break
        for fssp_key in dn_dlogdp_fssp_mean:
            if time == fssp_key.split("_")[-1]:
                buf[fssp_key] = dn_dlogdp_fssp_mean[fssp_key]
                break
        dn_dlogdp_comp_mean[time] = buf

    if plot_exact_dn_dlogdp is True:
        sam_bins = None
        for time in dn_dlogdp_comp:
            for key in dn_dlogdp_comp[time]:
                if "Droplet" in key:
                    sam_bins = sam_bins_droplet
                elif "Aerosol" in key:
                    sam_bins = sam_bins_aerosol

            if sam_bins is None:
                raise ValueError("ERROR: SAM bins not set properly")

            PacePlots.plot_pace_dn_dlogdp(dn_dlogdp_comp[time],
                                          sam_bins=sam_bins, cas_bins=cas_bins, fssp_bins=fssp_bins)
    if plot_mean_dn_dlogdp is True:
        sam_bins = None
        for time in dn_dlogdp_comp_mean:
            for key in dn_dlogdp_comp_mean[time]:
                if "Droplet" in key:
                    sam_bins = sam_bins_droplet
                elif "Aerosol" in key:
                    sam_bins = sam_bins_aerosol

            if sam_bins is None:
                raise ValueError("ERROR: SAM bins not set properly")

            PacePlots.plot_pace_dn_dlogdp(dn_dlogdp_comp_mean[time],
                                          sam_bins=sam_bins, cas_bins=cas_bins, fssp_bins=fssp_bins)

    pass
