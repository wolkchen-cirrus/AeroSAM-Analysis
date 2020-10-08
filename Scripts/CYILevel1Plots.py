from AirborneParticleAnalysis import level1to2
from AirborneParticleAnalysis import StandardLevel1Plots
from AirborneParticleAnalysis import StandardLevel2Plots
from scipy.signal import find_peaks
import numpy as np
from os import listdir
from matplotlib import pyplot as plt
from AirborneParticleAnalysis import common


ucass_number = 1            # 1 or 2
profile = "Up"              # Up or Down
dn_slices = []              # if empty then auto selection based on num. conc. peaks, or specify altitudes in [] (m)
conc_type = "Number"        # Exact as passed into function, see readme
strat_sizes = [0, 5, 10]

# Level 1 Plots
plot_conc_profile = 0       # Plot concentration
plot_dn_slices = 0          # Plot dn/dlog(Dp) slices
plot_strat_size = 0         # Plot Stratified Size

# Level 2 Plots
plot_rmar = 0
plot_cint_dn = 1

debug_plots = 0             # Debugging plots, only enable if admin
save_plots = 0
show_plots = 1

if __name__ == "__main__":
    print "============ Executing the Plotting Script for CYI SUA Level 1 Data ============\n"

    # Starting the import of the level 1 data, the script "convert_all_0to1.py" must be run first to obtain this level 1
    # data in the first instance.
    data_dir = "C:\\Users\\JGirdwood\\University of Hertfordshire\\AeroSAM - Documents\\Data\\2020\\20-03-10\\level_1"
    data_files = listdir(data_dir)                                                      # Getting files list
    data_dict = {}
    for file_name in data_files:
        print "INFO: Importing file - %s" % file_name
        key_name = file_name.split(".")[0]                                              # Dict key name is file name
        file_path = data_dir + "\\" + file_name                                         # Find full file path
        data_dict[key_name] = level1to2.import_level1(file_path)                        # Importing, uses "cPickle"

    # Level 1 number concentration plots

    index = 0
    fig_dict = {}

    for data in data_dict:

        window = int(common.read_setting("conc_window_size"))

        if plot_conc_profile:
            f, t = StandardLevel1Plots.level1_conc_plot(data_dict[data], conc_type="Mass", ucass_number=ucass_number)
            t = t.replace(" ", "_").replace("/", "").replace("\n", "_")\
                .replace(":", "").replace(")", "").replace("(", "")
            fig_dict[t] = f

        if profile == "Up":
            mask = data_dict[data].up_profile_mask
        elif profile == "Down":
            mask = data_dict[data].down_profile_mask
        else:
            raise ValueError("ERROR: Profile is either \"Up\" or \"Down\" as str")
        if ucass_number == 1:
            nc = data_dict[data].number_concentration1[np.where(mask[:, 0] == 1)]
            bin_centres = data_dict[data].bin_centres_dp_um1
        elif ucass_number == 2:
            nc = data_dict[data].number_concentration2[np.where(mask[:, 0] == 1)]
            bin_centres = data_dict[data].bin_centres_dp_um2
        else:
            raise ValueError("ERROR: UCASS number is either 1 or 2 passed as int")

        alt = data_dict[data].alt[np.where(mask[:, 0] == 1)]
        nc = nc / 1e6

        if not dn_slices:
            ncp_index, props = find_peaks(np.squeeze(nc), prominence=1)
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
            if not isinstance(dn_slices[0], list):
                raise ValueError("ERROR: Input dn_slices must be list of lists i.e. [[1,2,3],[1,2,3]]")
            alt_list = []
            try:
                dn = dn_slices[index]
            except IndexError:
                raise ValueError("ERROR: Number of input lists must be equal to number of flights")
            for dn_slice in dn:
                alt_list.append(level1to2.fetch_row(altitude=dn_slice, level1_data=data_dict[data], profile=profile))
        else:
            raise ValueError("ERROR: Invalid dn_slices variable, can only be list of floats or None")

        if plot_dn_slices:
            f, t = StandardLevel1Plots.level1_psd_plot(data_dict[data], list(alt_list), ucass_number=ucass_number)
            t = t.replace(" ", "_").replace("/", "").replace("\n", "_") \
                .replace(":", "").replace(")", "").replace("(", "")
            fig_dict[t] = f

        if plot_strat_size:
            f, t = StandardLevel1Plots.\
                level1_stratified_size(data_dict[data], strat_sizes, ucass_number=ucass_number, profile=profile)

        if plot_rmar:
            if data_dict[data].ucass_gain1 != data_dict[data].ucass_gain2:
                raise RuntimeError("ERROR: UCASS gains must be the same for plot")

            mask = np.add(data_dict[data].up_profile_mask, data_dict[data].down_profile_mask)

            dn1 = np.asarray(data_dict[data].dn_dlogdp1.values())[np.where(mask[:, 0] == 1), :]
            dn2 = np.asarray(data_dict[data].dn_dlogdp2.values())[np.where(mask[:, 0] == 1), :]

            dt = data_dict[data].datetime
            title_datetime = \
                dt[0:4] + "/" + dt[4:6] + "/" + dt[6:8] + " " + dt[8:10] + ":" + dt[10:12] + ":" + dt[12:14]

            dn1 = dn1.flatten()
            dn2 = dn2.flatten()
            new_dn1 = []
            new_dn2 = []
            for d1, d2 in zip(dn1, dn2):
                if d1 == d2 == 0:
                    pass
                else:
                    new_dn1.append(d1)
                    new_dn2.append(d2)
            dn1 = new_dn1
            dn2 = new_dn2

            rmar = level1to2.rma_regression(dn1, dn2)
            f, t = StandardLevel2Plots.plot_rebin_1to1(dn1, dn2, bin_centres, rmar,
                                                       ["UCASS 1", "UCASS 2", title_datetime], contour=None)
            t = t.replace(" ", "_").replace(")", "").replace("(", "").replace("/", "").replace(":", "")
            fig_dict[t] = f

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
