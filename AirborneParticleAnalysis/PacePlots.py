from matplotlib import pyplot as plt
from matplotlib import font_manager as plt_fnt
from matplotlib import rcParams as mplParams
from matplotlib import lines
from matplotlib.legend import Legend
import matplotlib.ticker as plticker
from AirborneParticleAnalysis import common


plt.style.use("ggplot")
prop = plt_fnt.FontProperties(family=['serif'])
mplParams["font.family"] = prop.get_name()
mplParams['hatch.linewidth'] = 0.5
mplParams['mathtext.default'] = "regular"


def plot_rebin_1to1(data_ref, data_sam, regression_data, mode):

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(8.3, 8.3))
    ax = fig.add_axes([0.2, 0.15, 0.75, 0.7])
    title_string = "Integrated dn/dlog(Dp), SUA vs CAS"
    ax.set_title(title_string, fontsize="small")
    ax.set_ylabel('CAS', fontsize="small")
    ax.set_xlabel('SUA in %s Mode' % mode, fontsize="small")

    x12 = regression_data[0]
    y12 = regression_data[1]
    r2 = regression_data[2]
    m = regression_data[5]

    marker_style = dict(linestyle='none', marker='x', markersize=5, fillstyle='none', color=(0, 0, 0), linewidth=0.7)

    ax.plot(data_sam, data_ref, **marker_style)

    marker_style["linestyle"] = '-.'
    marker_style["marker"] = ''
    ax.plot(x12, y12, **marker_style)
    h1 = lines.Line2D([], [], **marker_style)

    marker_style["linestyle"] = 'solid'
    ax.plot([0, x12[-1]], [0, x12[-1]], **marker_style)
    h2 = lines.Line2D([], [], **marker_style)

    if mode == "Droplet":
        ax.text(100, 1400, "$r^{2} = $%f" % r2, fontsize="small")
        ax.text(100, 1250, "$m = $%f" % m, fontsize="small")
    elif mode == "Aerosol":
        ax.text(7, 128, "$r^{2} = $%f" % r2, fontsize="small")
        ax.text(7, 114, "$m = $%f" % m, fontsize="small")

    ax.set_ylim(ymin=0)
    ax.set_xlim(xmin=0)

    if mode == "Droplet":
        loc = plticker.MultipleLocator(base=500)
    elif mode == "Aerosol":
        loc = plticker.MultipleLocator(base=50)
    else:
        raise ValueError("ERROR: Unrecognised mode %s" % mode)
    ax.yaxis.set_major_locator(loc)
    ax.xaxis.set_major_locator(loc)

    leg = Legend(ax, [h1, h2], ["Regression Line", "y = x"], frameon=False, fontsize="small", loc=2)
    ax.add_artist(leg)

    plt.show()

    return


def plot_pace_dn_dlogdp(data_dict, sam_bins=None, cas_bins=None, fssp_bins=None):

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(12, 8))
    ax = fig.add_axes([0.15, 0.2, 0.8, 0.7])
    time = common.seconds_to_timestamp(int(data_dict.keys()[0].split("_")[-1]))
    title_string = "PSD for the UCASS - V2 and Reference Instrumentation at %s" % time
    ax.set_title(title_string, fontsize="small")

    marker_styles = ["x", "o", "+", "s", "D"]
    line_style = ['solid', 'dotted', 'dashed', 'dashdot']
    marker_style = dict(linestyle=':', marker='x', markersize=3, fillstyle='none', color=(0, 0, 0), linewidth=0.7,
                        capsize=2)
    legend1_style = dict(marker='x', color=(0, 0, 0), linestyle='None', fillstyle='none', linewidth=0.7)

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
            legend1_style['marker'] = marker_styles[index]
            legend1_style['linestyle'] = line_style[index]
            patch1_handle = lines.Line2D([], [], **legend1_style)
            patch1_handles.append(patch1_handle)

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

            legend1_style['marker'] = marker_styles[index]
            legend1_style['linestyle'] = line_style[index]
            patch1_handle = lines.Line2D([], [], **legend1_style)
            patch1_handles.append(patch1_handle)

            marker_style['marker'] = marker_styles[index]
            marker_style['linestyle'] = line_style[index]
            line_handle = ax.errorbar(bins[cas_noise:], data[cas_noise:], yerr=dn_err[cas_noise:], **marker_style)
            line_handles.append(line_handle)

            leg_labels.append(leg_label)

            index += 1

    if bins is None:
        raise ValueError("ERROR: No bins specified")

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
