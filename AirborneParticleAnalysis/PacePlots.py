from matplotlib import pyplot as plt
from matplotlib import font_manager as plt_fnt
from matplotlib import rcParams as mplParams
from matplotlib import lines
from matplotlib.legend import Legend
import matplotlib.colors as mpl_col
import matplotlib.ticker as plticker
from matplotlib.ticker import FixedLocator, FixedFormatter
from AirborneParticleAnalysis import common
from AirborneParticleAnalysis import level1to2
import matplotlib as mpl
import numpy as np
import matplotlib.gridspec as gridspec
from string import ascii_lowercase
from scipy import stats


prop = plt_fnt.FontProperties(family=['serif'])
mplParams["font.family"] = prop.get_name()
mplParams['hatch.linewidth'] = 0.5
mplParams['mathtext.default'] = "regular"


def conc_pry_2020(conc_diff, pry):

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(8.3, 5.5))
    ax = fig.add_axes([0.2, 0.3, 0.75, 0.6])

    m, c, _, p, _ = stats.linregress(pry, conc_diff)

    ax.plot(pry, conc_diff, linestyle='none', marker='x', color=(0, 0, 0))
    ax.plot([min(pry), max(pry)], [m*min(pry)+c, m*max(pry)+c], color=(0, 0, 0), linestyle=':')
    ax.text(0.15, 0.9, 'p = %s' % round(p, 3), ha='center', va='center', transform=ax.transAxes, fontsize="small")

    ax.set_title("Concentration Error Against PRY Variation", fontsize="small")
    ax.set_ylabel("Number Concentration\nDifference " + r'($cm^{-3}$)', fontsize="small", labelpad=0.3)
    ax.set_xlabel("PRY Variation (AU)", fontsize="small")
    ax.set_ylim(0, 150)
    ax.yaxis.grid(True)

    ax.tick_params(axis="x", labelsize="x-small")
    ax.tick_params(axis="y", labelsize="x-small")

    return


def conc_asp_2020(conc_diff, asp):

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(8.3, 5.5))
    ax = fig.add_axes([0.2, 0.3, 0.75, 0.6])

    m, c, _, p, _ = stats.linregress(asp, conc_diff)

    ax.plot(asp, conc_diff, linestyle='none', marker='x', color=(0, 0, 0))
    ax.plot([min(asp), max(asp)], [m*min(asp)+c, m*max(asp)+c], color=(0, 0, 0), linestyle=':')
    ax.text(0.15, 0.9, 'p = %s' % round(p, 3), ha='center', va='center', transform=ax.transAxes, fontsize="small")

    ax.set_title("Concentration Error Against Airspeed", fontsize="small")
    ax.set_ylabel("Number Concentration\nDifference " + r'($cm^{-3}$)', fontsize="small", labelpad=0.3)
    ax.set_xlabel(r'Airspeed ($ms^{-1}$)', fontsize="small")
    ax.set_ylim(0, 150)
    ax.yaxis.grid(True)

    ax.tick_params(axis="x", labelsize="x-small")
    ax.tick_params(axis="y", labelsize="x-small")

    return


def eff_dia_plot_2020(eff_dia_diff, dt):

    k = 0
    for i, j, k in zip(dt[:-1], dt[1:], range(len(dt)-1)):
        num1 = float(i.split(":")[0] + i.split(":")[1])
        num2 = float(j.split(":")[0] + j.split(":")[1])
        if num2 < num1:
            break
        else:
            pass

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(8.3, 5.5))

    ax = fig.add_axes([0.2, 0.3, 0.75, 0.6])

    ax.plot(range(len(eff_dia_diff)), eff_dia_diff, linestyle='none', marker='x', color=(0, 0, 0))

    x_formatter = FixedFormatter(dt)
    x_locator = FixedLocator(range(len(eff_dia_diff)))
    ax.xaxis.set_major_formatter(x_formatter)
    ax.xaxis.set_major_locator(x_locator)
    ax.tick_params(axis="x", labelsize="x-small", labelrotation=60)
    ax.tick_params(axis="y", labelsize="x-small")

    ax.plot([-0.25, k+0.25], [5, 5], color='tab:blue')
    ax.plot([k+0.75, len(dt)-0.75], [5, 5], color='tab:orange')
    ax.plot([-0.25, -0.25], [5, 4.5], color='tab:blue')
    ax.plot([k+0.25, k+0.25], [5, 4.5], color='tab:blue')
    ax.plot([k + 0.75, k + 0.75], [5, 4.5], color='tab:orange')
    ax.plot([len(dt) - 0.75, len(dt) - 0.75], [5, 4.5], color='tab:orange')

    ax.text((-0.25 + k+0.25) / 2, 4.5, "2021-09-28", size='x-small', ha='center', va='top', color='tab:blue')
    ax.text((k+0.75 + len(dt)-0.75) / 2, 4.5, "2021-09-29", size='x-small', ha='center', va='top', color='tab:orange')

    ax.set_title("Effective Diameter Variation Over Campaign", fontsize="small")
    ax.set_ylabel("Effective Diameter\n(Talon - Static)", fontsize="small", labelpad=0.3)
    ax.set_xlabel("Take-Off Time (hh:mm)", fontsize="small")
    ax.set_ylim(-6, 6)
    ax.yaxis.grid(True)

    return


def simple_correlation(data_y, data_x, label_y, label_x, title, regress=False):
    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(12, 30))
    ax = fig.add_axes()
    ax.set_title(title, fontsize="small")
    ax.set_ylabel(label_y, fontsize="small")
    ax.set_xlabel(label_x, fontsize="small")
    ax.scatter(data_y, data_x, s=0.3)

    if regress is True:
        [x12, y12, r2, _, _, m] = level1to2.rma_regression(data_x, data_y)

        ax.text(plt.xlim()[0], plt.ylim()[1], r2)
        ax.text(plt.xlim()[0], plt.ylim()[1] - 100, m)

        ax.plot(x12, y12, linestyle='-.', marker='', color=(0, 0, 0))
        ax.plot([0, x12[-1]], [0, x12[-1]], linestyle='solid', marker='', color=(0, 0, 0))
    return


def plot_pace_dn_dlogdp_2020(talon_data, static_data, talon_bins, static_bins,
                             dt_string_arr, d_eff_s, d_eff_t):

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(12, 18))
    fig.suptitle('Data for Flights on %s' % dt_string_arr[0].split(" ")[0])

    col_plots = 2
    row_plots = int(np.ceil(len(talon_data)/float(col_plots)))
    gs = gridspec.GridSpec(ncols=col_plots, nrows=row_plots+1, figure=fig)

    ax_dict = {}
    for t, s, index, dt, de_s, de_t in \
            zip(talon_data, static_data, range(len(talon_data)), dt_string_arr, d_eff_s, d_eff_t):

        ax_dict[index] = fig.add_subplot(gs[index / col_plots, index % col_plots], xscale='log')

        talon_widths = [-1*(j-i) for i, j in zip(talon_bins[:-1], talon_bins[1:])]
        static_widths = [-1*(j-i) for i, j in zip(static_bins[:-1], static_bins[1:])]

        ax_dict[index].bar(talon_bins[1:], t[0], width=talon_widths, align="edge", alpha=0.5,
                           edgecolor="none", color='tab:blue')
        ax_dict[index].bar(static_bins[1:], s[0], width=static_widths, align="edge", alpha=0.5,
                           edgecolor="none", color='tab:orange')

        y_max = 2000
        ax_dict[index].plot([de_s, de_s], [0, y_max], linestyle='-.', marker='', color=(1, 0, 0))
        ax_dict[index].plot([de_t, de_t], [0, y_max], linestyle='-.', marker='', color=(0, 0, 1))

        title_string = str(dt).replace(' ', '\n')
        ax_dict[index].text(1, 150, title_string, fontsize='x-small')
        ax_dict[index].text(1, 1700, ascii_lowercase[index] + ")", fontsize='small')
        ax_dict[index].set_ylim(ymin=0, ymax=y_max)

        for tick in ax_dict[index].xaxis.get_major_ticks():
            tick.label.set_fontsize('x-small')
        for tick in ax_dict[index].yaxis.get_major_ticks():
            tick.label.set_fontsize('x-small')
        ax_dict[index].set_yticklabels(ax_dict[index].get_yticks(), rotation=90)
        ax_dict[index].yaxis.set_major_formatter(plticker.FormatStrFormatter('%5d'))
        ax_dict[index].tick_params(axis='both', which='major', pad=0.5)

    ax_dict[ax_dict.keys()[0]].set_title("Upwards Profile", fontsize="small")
    ax_dict[ax_dict.keys()[1]].set_title("Downwards Profile", fontsize="small")

    ax_dict[ax_dict.keys()[-1]].set_xlabel(r'Particle Diameter ($\mu m$)', fontsize="x-small", labelpad=0.3)
    ax_dict[ax_dict.keys()[-2]].set_xlabel(r'Particle Diameter ($\mu m$)', fontsize="x-small", labelpad=0.3)
    for i in range(len(talon_data)):
        if i % 2 is 1:
            pass
        else:
            ax_dict[ax_dict.keys()[i]].set_ylabel(r'dN/dlog($D_{p}$) ($cm^{3}$)', fontsize="x-small", labelpad=0.5)

    ax_dict[row_plots*2] = fig.add_subplot(gs[-1, :])
    ax_dict[row_plots*2].axis('off')
    leg_handles = [ax_dict[row_plots*2].fill_between([], [], [], color='tab:blue', alpha=0.5, edgecolor="none"),
                   ax_dict[row_plots*2].fill_between([], [], [], color='tab:orange', alpha=0.5, edgecolor="none"),
                   lines.Line2D([], [], linestyle='-.', marker='', color=(0, 0, 1)),
                   lines.Line2D([], [], linestyle='-.', marker='', color=(1, 0, 0))]
    leg_names = ['Talon UCASS', 'Static UCASS', 'Effective Diameter Talon', 'Effective Diameter Static']
    leg = Legend(ax_dict[row_plots*2], leg_handles, leg_names, frameon=False, fontsize="small", ncol=2, loc=8)
    ax_dict[row_plots*2].add_artist(leg)

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


def _add_interval(ax, x_data, y_data, caps="  "):
    line = ax.add_line(lines.Line2D(x_data, y_data))
    anno_args = {
        'ha': 'center',
        'va': 'center',
        'color': line.get_color()
    }
    a0 = ax.annotate(caps[0], xy=(x_data[0], y_data[0]), **anno_args)
    a1 = ax.annotate(caps[1], xy=(x_data[1], y_data[1]), **anno_args)
    return line, (a0, a1)
