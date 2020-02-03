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
from os import listdir


# Booleans for specifying which plots are required in the analysis
plot_mean_dn_dlogdp = True      # dn/dlog(Dp) averaged over a height/time specified in settings.txt
plot_exact_dn_dlogdp = False    # dn/dlog(Dp) at an exact point with no average.
plot_rebin_1to1 = True          # Re-binned data 1 to 1 plot.


if __name__ == "__main__":
    print "============ Executing the Plotting Script for the PaCE 2019 Campaign Data ============\n"

    # Starting the import of the level 1 data, the script "convert_all_0to1.py" must be run first to obtain this level 1
    # data in the first instance.
    data_dir = "C:\\Users\\jg17acv\\University of Hertfordshire\\AeroSAM - Documents\\Data\\2019\\19-09-28\\level_1"
    station_altitude_asl_mm = float(common.read_setting("station_altitude_asl_mm"))     # Altitude of CAS above SL
    data_files = listdir(data_dir)                                                      # Getting files list
    data_dict = {}
    for file_name in data_files:
        print "INFO: Importing file - %s" % file_name
        key_name = file_name.split(".")[0]                                              # Dict key name is file name
        file_path = data_dir + "\\" + file_name                                         # Find full file path
        data_dict[key_name] = level1to2.import_level1(file_path)                        # Importing, uses "cPickle"

    # Pre-allocating space for the variables used in analysis, this improves the running speed of the script.
    fssp_bins = None                # FSSP bin centres
    cas_bins = None                 # CAS bin centres
    sam_bins_droplet = None         # UCASS bin centres (droplet mode, low gain)
    sam_bins_aerosol = None         # UCASS bin centres (Aerosol mode, high gain)
    sam_bbs_droplet = None          # UCASS bin boundaries (droplet mode, low gain)
    sam_bbs_aerosol = None          # UCASS bin boundaries (Aerosol mode, high gain)
    dn_dlogdp_sam = {}              # dn/dlog(Dp) measured by UCASS at specific altitude
    dn_dlogdp_cas = {}              # dn/dlog(Dp) measured by CAS at specific time
    dn_dlogdp_fssp = {}             # dn/dlog(Dp) measured by FSSP at specific time
    dn_dlogdp_sam_mean = {}         # dn/dlog(Dp) measured by UCASS at mean altitude
    dn_dlogdp_cas_mean = {}         # dn/dlog(Dp) measured by CAS at mean time
    dn_dlogdp_fssp_mean = {}        # dn/dlog(Dp) measured by FSSP at mean time
    dn_rebin_sam = {}               # Integrated, discretised dn/dlog(Dp) for UCASS in it's own gain mode
    dn_rebin_cas_droplet = {}       # Integrated, discretised dn/dlog(Dp) for CAS according to UCASS droplet bins
    dn_rebin_fssp_droplet = {}      # Integrated, discretised dn/dlog(Dp) for FSSP according to UCASS droplet bins
    dn_rebin_cas_aerosol = {}       # Integrated, discretised dn/dlog(Dp) for CAS according to UCASS aerosol bins
    dn_rebin_fssp_aerosol = {}      # Integrated, discretised dn/dlog(Dp) for FSSP according to UCASS aerosol bins
    times = []                      # List of measurement times used to select appropriate CAS/FSSP data

    # The first part of the analysis determines the correct altitude, or range of altitudes, the SAM data is taken from.
    # It also fills up the relevant data buffers, and calculates the time in UTC that the data was recorded. This is
    # then used in the following analysis pass for CAS and FSSP. Note the SAM data is in UTC+1h and the CAS data is in
    # UTC+0h, so 3600 seconds needed to be added to CAS time.
    print "INFO: Analysis pipeline first pass"
    # Cycle through all data
    for key in data_dict:
        if ("CAS" in key) or ("FSSP" in key):
            # Skip the time based data on the first round, this loop calculates the times of measurement.
            print "INFO: Skipping CAS data and FSSP data on first round"
            continue
        else:
            print "INFO: Processing SAM data"

            # Find dn/dlog(Dp), mean dn/dlog(Dp), and the matching altitude (stored in "row_number")
            dn_buf = data_dict[key].dn_dlogdp
            row_number = level1to2.fetch_row(altitude=station_altitude_asl_mm, level1_data=data_dict[key])[0]
            rows_for_mean = level1to2.fetch_row_tolerance(altitude=station_altitude_asl_mm, level1_data=data_dict[key])
            mean_dn_dlogdp = level1to2.mean_dn_dlogdp(data_dict[key], rows_for_mean[0])

            # Get the matching time
            current_time = level1to2.get_time_from_alt(data_dict[key], row_number) + 3600   # 1hour time difference
            times.append(current_time)

            # Assigning the bin centres and boundaries for droplet or aerosol mode. Also re-bins dn/dlog(Dp) for a
            # 1-to-1 size dependant comparison.
            tags = data_dict[key].tags
            if "Droplet" in tags:
                sam_bins_droplet = data_dict[key].bin_centres_dp_um
                sam_bbs_droplet = data_dict[key].bin_bounds_dp_um
                dn_rebin = level1to2.rebin_dn_dlogdp(mean_dn_dlogdp[0], sam_bins_droplet, sam_bbs_droplet[:10])
                dn_key = key + "_" + "Droplet"      # Assigning this flag to the key makes it easy to match up later
            elif "Aerosol" in tags:
                sam_bins_aerosol = data_dict[key].bin_centres_dp_um
                sam_bbs_aerosol = data_dict[key].bin_bounds_dp_um
                dn_rebin = level1to2.rebin_dn_dlogdp(mean_dn_dlogdp[0], sam_bins_aerosol, sam_bbs_aerosol)
                dn_key = key + "_" + "Aerosol"      # Assigning this flag to the key makes it easy to match up later
            else:
                raise AttributeError("ERROR: Gain mode not specified")

            # Assigning data to pre-allocated storage buffers
            dn_key = dn_key + "_" + str(current_time)
            dn_dlogdp_sam[dn_key] = dn_buf[row_number]
            dn_dlogdp_sam_mean[dn_key] = mean_dn_dlogdp
            dn_rebin_sam[dn_key] = dn_rebin

    # The second pass is very similar to the first, except it calculates quantities for the time based data (e.g. the
    # FSSP and CAS data) using the flight times calculated in the first loop.
    print "INFO: Analysis pipeline Second pass"
    for time in times:              # Looping through times, since a different value needs to be matched each SAM datum
        for key in data_dict:       # For each time, loop through all level 1 data
            if ("CAS" in key) or ("FSSP" in key):

                # Running code which is common to both instruments, level 2 quantities are calculated in the same way
                # for all time based data.
                print "INFO: Processing CAS and FSSP data"
                dn_buf = data_dict[key].dn_dlogdp
                row_number = level1to2.fetch_row(time=time, level1_data=data_dict[key])[0]
                rows_for_mean = level1to2.fetch_row_tolerance(time=time, level1_data=data_dict[key])
                mean_dn_dlogdp = level1to2.mean_dn_dlogdp(data_dict[key], rows_for_mean[0])
                dn_key = key + "_" + str(time)

                # Assigning CAS specific variables
                if "CAS" in key:
                    cas_bins = data_dict[key].bin_centres_dp_um
                    dn_dlogdp_cas[dn_key] = dn_buf[row_number]
                    dn_dlogdp_cas_mean[dn_key] = mean_dn_dlogdp
                    dn_rebin_cas_aerosol[dn_key] = level1to2.rebin_dn_dlogdp(mean_dn_dlogdp[0], cas_bins,
                                                                             sam_bbs_aerosol)
                    dn_rebin_cas_droplet[dn_key] = level1to2.rebin_dn_dlogdp(mean_dn_dlogdp[0], cas_bins,
                                                                             sam_bbs_droplet[:10])

                # Assigning FSSP specific variables
                elif "FSSP" in key:
                    fssp_bins = data_dict[key].bin_centres_dp_um
                    dn_dlogdp_fssp[dn_key] = dn_buf[row_number]
                    dn_dlogdp_fssp_mean[dn_key] = mean_dn_dlogdp
                    dn_rebin_fssp_aerosol[dn_key] = level1to2.rebin_dn_dlogdp(mean_dn_dlogdp[0], fssp_bins,
                                                                              sam_bbs_aerosol)
                    dn_rebin_fssp_droplet[dn_key] = level1to2.rebin_dn_dlogdp(mean_dn_dlogdp[0], fssp_bins,
                                                                              sam_bbs_droplet[:10])
            else:
                # Skip SAM data because this has already been processed at this point in the script.
                print "INFO: Skipping SAM data on second round"
                continue

    # This part of the script matches the appropriate CAS/FSSP and SAM data together in a dictionary for easy plotting
    print "INFO: Analysis pipeline third pass"
    dn_dlogdp_comp = {}
    dn_dlogdp_comp_mean = {}
    # Cycling through the SAM data (iterative search algorithm)
    for sam_key in dn_dlogdp_sam:
        time = sam_key.split("_")[-1]               # Get time in seconds for matching with CAS/FSSP
        buf = {sam_key: dn_dlogdp_sam[sam_key]}     # Assign SAM data to buffer under same key
        # Loop through comparative data for each datum in the dn_dlogdp_sam dict
        for cas_key in dn_dlogdp_cas:
            if time == cas_key.split("_")[-1]:          # Find which time matches
                buf[cas_key] = dn_dlogdp_cas[cas_key]   # Add to buffer
                break                                   # If a match is found, stop searching
        for fssp_key in dn_dlogdp_fssp:                 # Repeat with FSSP
            if time == fssp_key.split("_")[-1]:
                buf[fssp_key] = dn_dlogdp_fssp[fssp_key]
                break
        dn_dlogdp_comp[time] = buf                      # Assign buf to storage dictionary for plotting

    # Exactly the same process as above but for the mean dn/dlog(Dp) data
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

    # Finding comparative data for the dn/dlog(Dp) re-binned data. I.e. this is the process of matching the correct SAM
    # data to the correct FSSP and CAS data, re-binned into the correct UCASS gain mode boundaries.
    drop_buf_cas = []
    drop_buf_sam = []
    aerosol_buf_sam = []
    aerosol_buf_cas = []
    for sam_key in dn_rebin_sam:
        time = sam_key.split("_")[-1]
        for cas_key in dn_rebin_cas_aerosol:
            if time == cas_key.split("_")[-1]:
                if "Aerosol" in sam_key:
                    aerosol_buf_sam.append(dn_rebin_sam[sam_key])
                    aerosol_buf_cas.append(dn_rebin_cas_aerosol[cas_key])
        for cas_key in dn_rebin_cas_droplet:
            if time == cas_key.split("_")[-1]:
                if "Droplet" in sam_key:
                    drop_buf_sam.append(dn_rebin_sam[sam_key])
                    drop_buf_cas.append(dn_rebin_cas_droplet[cas_key])
    # The above gives us lists of lists so we need to flatten them (below).
    drop_buf_cas = [val for sublist in drop_buf_cas for val in sublist]
    drop_buf_sam = [val for sublist in drop_buf_sam for val in sublist]
    aerosol_buf_cas = [val for sublist in aerosol_buf_cas for val in sublist]
    aerosol_buf_sam = [val for sublist in aerosol_buf_sam for val in sublist]
    # Getting linear regression line and values.
    rebin_drop_regression = level1to2.rma_regression(drop_buf_sam, drop_buf_cas)
    rebin_aerosol_regression = level1to2.rma_regression(aerosol_buf_sam, aerosol_buf_cas)

    # Making the 1 to 1 plot of the re-binned data
    if plot_rebin_1to1 is True:
        PacePlots.plot_rebin_1to1(drop_buf_cas, drop_buf_sam, rebin_drop_regression, "Droplet")
        PacePlots.plot_rebin_1to1(aerosol_buf_cas, aerosol_buf_sam, rebin_aerosol_regression, "Aerosol")

    # Plotting the dn/dlog(Dp) graphs for the exact times and altitudes
    if plot_exact_dn_dlogdp is True:
        sam_bins = None
        for time in dn_dlogdp_comp:
            for key in dn_dlogdp_comp[time]:        # Search keys to find which gain bin centres to plot against
                if "Droplet" in key:
                    sam_bins = sam_bins_droplet     # If droplet gain mode
                elif "Aerosol" in key:
                    sam_bins = sam_bins_aerosol     # If aerosol gain mode

            if sam_bins is None:
                raise ValueError("ERROR: SAM bins not set properly")

            PacePlots.plot_pace_dn_dlogdp(dn_dlogdp_comp[time],
                                          sam_bins=sam_bins, cas_bins=cas_bins, fssp_bins=fssp_bins)

    # Exactly the same as above but for the mean dn/dlog(Dp)
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
