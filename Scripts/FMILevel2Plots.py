from AirborneParticleAnalysis import level1to2
from AirborneParticleAnalysis import StandardLevel1Plots
from AirborneParticleAnalysis import StandardLevel2Plots
from scipy.signal import find_peaks
import numpy as np
from os import listdir
from matplotlib import pyplot as plt
from AirborneParticleAnalysis import common


profile = "Up"              # Up or Down
dn_slices = [400, 500, 600, 700, 800, 900]          # if empty then auto selection based on num. conc. peaks, or specify altitudes in [] (m)
conc_type = "Number"        # Exact as passed into function, see readme
strat_sizes = [0, 5, 10]

# Level 1 Plots
plot_conc_profile = 1       # Plot concentration
plot_dn_slices = 0          # Plot dn/dlog(Dp) slices
plot_strat_size = 0         # Plot Stratified Size

# Level 2 Plots
plot_rmar = 0
plot_cint_dn = 0

debug_plots = 0             # Debugging plots, only enable if admin
save_plots = 0
show_plots = 1

if __name__ == "__main__":
    print "============ Executing the Plotting Script for FMI SUA Level 2 Data ============\n"

    # Starting the import of the level 2 data, the script "convert_all_0to1.py" must be run first to obtain this level 1
    # data in the first instance.
    data_dir = common.read_setting("base_data_dir")
    t1 = listdir(data_dir)
    data_files = []
    for date in t1:
        date_dir = data_dir + "\\" + date + "\\" + "level_2"
        t2 = listdir(date_dir)
        for fnm in t2:
            if "FMITalon_" in fnm:
                data_files.append(date_dir + "\\" + fnm)

    data_dict = {}
    for file_name in data_files:
        print "INFO: Importing file - %s" % file_name
        key_name = file_name.split("\\")[-1].split(".")[0]                              # Dict key name is file name
        data_dict[key_name] = level1to2.import_level1(file_name)                        # Importing, uses "cPickle"

    # Level 1 number concentration plots

    index = 0
    fig_dict = {}

    for data in data_dict:

        window = int(common.read_setting("conc_window_size"))

        if plot_conc_profile:

            profile_number = data_dict[data].up_profile_mask.shape[1]
            for prof in range(profile_number):
                f, t = StandardLevel2Plots.level2_conc_plot(data_dict[data],
                                                            prof_num=prof, conc_type=conc_type, dn=dn_slices,
                                                            asp_lim=(0, 30), conc_lim=(0, 1500), dn_lim=(0, 450),
                                                            alt_lim=(0, 1350))
                t = t.replace(" ", "_").replace("/", "").replace("\n", "_")\
                    .replace(":", "").replace(")", "").replace("(", "")
                fig_dict[t] = f

        if profile == "Up":
            mask = data_dict[data].up_profile_mask
        elif profile == "Down":
            mask = data_dict[data].down_profile_mask
        else:
            raise ValueError("ERROR: Profile is either \"Up\" or \"Down\" as str")
        nc = data_dict[data].number_concentration[np.where(mask[:, 0] == 1)]
        bin_centres = data_dict[data].bin_centres_dp_um

        alt = data_dict[data].alt[np.where(mask[:, 0] == 1)]
        nc = nc / 1e6

        if not dn_slices:
            ncp_index, props = find_peaks(np.squeeze(nc), prominence=5)
            proms = props["prominences"]
            tops = sorted(zip(proms, ncp_index), reverse=True)[:4]
            top_index = [j for i, j in tops]
            alt_list = alt[top_index]
            alt_list = [alt[0] for alt in alt_list]
            if debug_plots:
                plt.plot(nc)
                plt.plot(top_index, nc[top_index], 'x')
                pass
        elif dn_slices:
            alt_list = []
            for dn_slice in dn_slices:
                alt_list.append(level1to2.fetch_row(altitude=dn_slice, level1_data=data_dict[data], profile=profile))
        else:
            raise ValueError("ERROR: Invalid dn_slices variable, can only be list of floats or None")

        if plot_dn_slices:
            f, t = StandardLevel1Plots.level1_psd_plot(data_dict[data], list(alt_list))
            t = t.replace(" ", "_").replace("/", "").replace("\n", "_") \
                .replace(":", "").replace(")", "").replace("(", "")
            fig_dict[t] = f

        if plot_strat_size:
            f, t = StandardLevel1Plots.\
                level1_stratified_size(data_dict[data], strat_sizes, profile=profile)

        if plot_cint_dn:
            f, t = StandardLevel2Plots.plot_cint_dn(data_dict[data], end_swipe=100)
            t = t.replace(" ", "_").replace(")", "").replace("(", "").replace("/", "").replace(":", "")\
                .replace("\n", "_")
            fig_dict[t] = f

    index += 1

    if save_plots:
        for t in fig_dict:
            name = data_dir.replace("level_1", "") + t + ".pdf"
            fig_dict[t].savefig(name, format='pdf')

    if show_plots:
        plt.show()
