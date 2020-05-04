from AirborneParticleAnalysis import level1to2
from AirborneParticleAnalysis import StandardLevel1Plots
from scipy.signal import find_peaks
import numpy as np
from os import listdir
from matplotlib import pyplot as plt
from AirborneParticleAnalysis import common


ucass_number = 1            # 1 or 2
profile = "Up"              # Up or Down
dn_slices = []              # if empty then auto selection based on num. conc. peaks, or specify altitudes in [] (m)
conc_type = "Number"        # Exact as passed into function, see readme

plot_conc_profile = 0       # Plot concentration
plot_dn_slices = 0          # Plot dn/dlog(Dp) slices

debug_plots = 1             # Debugging plots, only enable if admin

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
    for data in data_dict:

        alt = data_dict[data].alt
        window = int(common.read_setting("conc_window_size"))

        if plot_conc_profile:
            StandardLevel1Plots.level1_conc_plot(data_dict[data], conc_type="Number", ucass_number=ucass_number)

        if ucass_number == 1:
            if profile == "Up":
                mask = data_dict[data].up_profile_mask
                nc = data_dict[data].number_concentration1[np.where(mask[:, 0] == 1)]
                nc = np.convolve(nc[:, 0], np.ones((window,)) / window, mode="same")
            elif profile == "Down":
                mask = data_dict[data].down_profile_mask
                nc = data_dict[data].number_concentration1[np.where(mask[:, 0] == 1)]
            else:
                raise ValueError("ERROR: Profile is either \"Up\" or \"Down\" as str")
        elif ucass_number == 2:
            if profile == "Up":
                mask = data_dict[data].up_profile_mask
                nc = data_dict[data].number_concentration2[np.where(mask[:, 0] == 1)]
            elif profile == "Down":
                mask = data_dict[data].down_profile_mask
                nc = data_dict[data].number_concentration2[np.where(mask[:, 0] == 1)]
            else:
                raise ValueError("ERROR: Profile is either \"Up\" or \"Down\" as str")
        else:
            raise ValueError("ERROR: UCASS number is either 1 or 2 passed as int")

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
            alt_list = []
            for dn_slice in dn_slices:
                alt_list.append(level1to2.fetch_row(altitude=dn_slice, level1_data=data_dict[data], profile=profile))
        else:
            raise ValueError("ERROR: Invalid dn_slices variable, can only be list of floats or None")

        if plot_dn_slices:
            StandardLevel1Plots.level1_psd_plot(data_dict[data], list(alt_list), ucass_number=ucass_number)

plt.show()
