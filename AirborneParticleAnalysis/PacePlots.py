from matplotlib import pyplot as plt
from matplotlib import font_manager as plt_fnt
from matplotlib import rcParams as mplParams
from matplotlib import lines
from matplotlib.legend import Legend
import matplotlib.colors as mpl_col
import matplotlib.ticker as plticker
from AirborneParticleAnalysis import common
import matplotlib as mpl
import numpy as np


# plt.style.use("ggplot")
prop = plt_fnt.FontProperties(family=['serif'])
mplParams["font.family"] = prop.get_name()
mplParams['hatch.linewidth'] = 0.5
mplParams['mathtext.default'] = "regular"


def plot_pace_dn_dlogdp_2020(talon_data, static_data, talon_bins, static_bins, dt_string_arr):

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(12, 12))

    col_plots = 2
    row_plots = int(np.ceil(len(talon_data)/float(col_plots)))

    ax_dict = {}
    for t, s, index, dt in zip(talon_data, static_data, range(len(talon_data)), dt_string_arr):

        ax_dict[index] = fig.add_subplot(row_plots, col_plots, index+1)

        talon_widths = [-1*(j-i) for i, j in zip(talon_bins[:-1], talon_bins[1:])]
        static_widths = [-1*(j-i) for i, j in zip(static_bins[:-1], static_bins[1:])]

        ax_dict[index].bar(talon_bins[1:], t[0], width=talon_widths, align="edge", alpha=0.5, label='Talon UCASS')
        ax_dict[index].bar(static_bins[1:], s[0], width=static_widths, align="edge", alpha=0.5, label='Static UCASS')

        title_string = "Take-off at %s" % dt
        ax_dict[index].set_title(title_string, fontsize="small")
        ax_dict[index].set_ylabel('Normalised Concentration\n(dN/dlogDp)', fontsize="small")
        ax_dict[index].set_xlabel(r'Particle Diameter ($\mu m$)', fontsize="small")
        ax_dict[index].set_ylim(ymin=0, ymax=2000)

        ax_dict[index].legend(frameon=False, fontsize="small")

        plt.xscale('log')

    plt.show()
    return


def plot_rebin_1to1(data_ref, data_sam, regression_data, mode, bin_centres):

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(8.3, 11.3))
    ax = fig.add_axes([0.2, 0.3, 0.75, 0.55])
    title_string = "Integrated dN/dlog(Dp), CAS vs SUA"
    ax.set_title(title_string, fontsize="small")
    ax.set_ylabel('CAS', fontsize="small")
    ax.set_xlabel('SUA in %s Mode' % mode, fontsize="small")

    x12 = regression_data[0]
    y12 = regression_data[1]
    r2 = regression_data[2]
    m = regression_data[5]

    data_z = []
    index = 0
    while True:
        if len(data_z) == len(data_sam):
            break
        data_z.append(bin_centres[index])
        index += 1
        if index == 16:
            index = 0

    marker_style = dict(linestyle='none', marker='x', markersize=5, fillstyle='none', color=(0, 0, 0), linewidth=0.7)
    cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = mpl_col.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = np.linspace(bin_centres[0], bin_centres[-1], len(bin_centres))
    norm = mpl_col.BoundaryNorm(bounds, cmap.N)

    ax.scatter(data_sam, data_ref, c=data_z, cmap=cmap, norm=norm)

    ax2 = fig.add_axes([0.2, 0.11, 0.75, 0.03])
    mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional',
                              ticks=np.linspace(bin_centres[0], bin_centres[-1], len(bin_centres)/2), boundaries=bounds,
                              format='%1i', orientation='horizontal')
    ax2.set_xlabel(r'Particle Diameter ($\mu m$)', fontsize="small")

    marker_style["linestyle"] = '-.'
    marker_style["marker"] = ''
    marker_style["color"] = 'black'

    if mode == "Droplet":
        # ax.text(15, 1000, "$r^{2} = $%f" % r2, fontsize="small")
        # ax.text(15, 900, "$m = $%f" % m, fontsize="small")

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_ylim(ymin=10)
        ax.set_xlim(xmin=10)

        reg_x = np.linspace(x12[0]+0.000001, x12[1], 1000)
        reg_y = np.linspace(y12[0], y12[1], 1000)

        ax.plot(reg_x, reg_y, **marker_style)
        h1 = lines.Line2D([], [], **marker_style)

    elif mode == "Aerosol":
        # ax.text(7, 128, "$r^{2} = $%f" % r2, fontsize="small")
        # ax.text(7, 114, "$m = $%f" % m, fontsize="small")

        ax.set_ylim(ymin=0)
        ax.set_xlim(xmin=0)

        ax.plot(x12, y12, **marker_style)
        h1 = lines.Line2D([], [], **marker_style)

    marker_style["linestyle"] = 'solid'
    ax.plot([0, x12[-1]], [0, x12[-1]], **marker_style)
    h2 = lines.Line2D([], [], **marker_style)

    if mode == "Droplet":
        pass
        # loc = plticker.MultipleLocator(base=10)
    elif mode == "Aerosol":
        loc = plticker.MultipleLocator(base=50)
        ax.yaxis.set_major_locator(loc)
        ax.xaxis.set_major_locator(loc)
    else:
        raise ValueError("ERROR: Unrecognised mode %s" % mode)

    legend1_style = dict(marker='None', linestyle='None', fillstyle='none')
    patches = [lines.Line2D([], [], **legend1_style), lines.Line2D([], [], **legend1_style)]

    leg = Legend(ax, [h1, h2, patches[0], patches[1]], ["Regression Line", "y = x", "$r^{2} = $%f" % r2, "$m = $%f" % m]
                 , frameon=False, fontsize="small", loc=4)
    ax.add_artist(leg)

    plt.show()

    return


def plot_pace_dn_dlogdp(data_dict, sam_bins=None, cas_bins=None, fssp_bins=None, plot_fssp=False):

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(12, 8))
    ax = fig.add_axes([0.15, 0.2, 0.8, 0.7])
    time = common.seconds_to_timestamp(int(data_dict.keys()[0].split("_")[-1]))
    title_string = "PSD for the UCASS - V2 and Reference Instrumentation at %s" % time
    ax.set_title(title_string, fontsize="small")

    marker_styles = ["x", "o", "+", "s", "D"]
    line_style = ['solid', 'dotted', 'dashed', 'dashdot']
    marker_style = dict(linestyle=':', marker='x', markersize=3, fillstyle='none', color='C0', linewidth=0.7,
                        capsize=2)
    legend1_style = dict(marker='x', color='C0', linestyle='None', fillstyle='none', linewidth=0.7)

    patch1_handles = []
    line_handles = []
    leg_labels = []
    index = 0
    bins = None
    cas_noise = int(common.read_setting("cas_noise_threshold_bin"))

    for key in data_dict:

        data = data_dict[key][0]
        dn_err = data_dict[key][1]

        if "SAM" in key:
            bins = sam_bins
            legend1_style['color'] = "C" + str(index)
            legend1_style['marker'] = marker_styles[index]
            legend1_style['linestyle'] = line_style[index]
            patch1_handle = lines.Line2D([], [], **legend1_style)
            patch1_handles.append(patch1_handle)

            marker_style['color'] = "C" + str(index)
            marker_style['marker'] = marker_styles[index]
            marker_style['linestyle'] = line_style[index]
            line_handle = ax.errorbar(bins, data, yerr=dn_err, **marker_style)
            line_handles.append(line_handle)

            leg_label = "SUA with UCASS-V2"
            leg_labels.append(leg_label)

            index += 1

    if bins is None:
        raise ValueError("ERROR: No bins specified")

    bins = None
    for key in data_dict:

        data = data_dict[key][0]
        dn_err = data_dict[key][1]

        if "CAS" in key:

            bins = cas_bins
            leg_label = "CAS at %sm" % str(int(common.read_setting("station_altitude_asl_mm"))/1000)

            legend1_style['color'] = "C" + str(index)
            legend1_style['marker'] = marker_styles[index]
            legend1_style['linestyle'] = line_style[index]
            patch1_handle = lines.Line2D([], [], **legend1_style)
            patch1_handles.append(patch1_handle)

            marker_style['color'] = "C" + str(index)
            marker_style['marker'] = marker_styles[index]
            marker_style['linestyle'] = line_style[index]
            line_handle = ax.errorbar(bins[cas_noise:], data[cas_noise:], yerr=dn_err[cas_noise:], **marker_style)
            line_handles.append(line_handle)

            leg_labels.append(leg_label)

            index += 1

    if bins is None:
        raise ValueError("ERROR: No bins specified")

    if plot_fssp:
        bins = None
        for key in data_dict:

            data = data_dict[key][0]
            dn_err = data_dict[key][1]

            if "FSSP" in key:

                bins = fssp_bins
                leg_label = "FSSP at %sm" % str(int(common.read_setting("station_altitude_asl_mm")) / 1000)

                legend1_style['marker'] = marker_styles[index]
                legend1_style['linestyle'] = line_style[index]
                patch1_handle = lines.Line2D([], [], **legend1_style)
                patch1_handles.append(patch1_handle)

                marker_style['marker'] = marker_styles[index]
                marker_style['linestyle'] = line_style[index]
                line_handle = ax.errorbar(bins, data, yerr=dn_err, **marker_style)
                line_handles.append(line_handle)

                leg_labels.append(leg_label)

                index += 1

    if bins is None:
        raise ValueError("ERROR: No bins specified")

    plt.xscale('log')

    ax.set_ylabel('Normalised Concentration\n(dN/dlogDp)', fontsize="small")
    ax.set_xlabel(r'Particle Diameter ($\mu m$)', fontsize="small")
    ax.set_ylim(ymin=0)
    ax.set_ylim(ymax=650)
    ax.set_xlim(xmin=0.5)

    leg = Legend(ax, patch1_handles, leg_labels, frameon=False, fontsize="small")
    ax.add_artist(leg)

    plt.show()
    return
