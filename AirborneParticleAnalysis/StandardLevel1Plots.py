from matplotlib import pyplot as plt
from matplotlib import font_manager as plt_fnt
from matplotlib import rcParams as mplParams
from matplotlib import lines
from matplotlib.legend import Legend
# import matplotlib.ticker as plticker
from AirborneParticleAnalysis import common
from AirborneParticleAnalysis import level1to2
import numpy as np


# plt.style.use("ggplot")
prop = plt_fnt.FontProperties(family=['serif'])
mplParams["font.family"] = prop.get_name()
mplParams['hatch.linewidth'] = 0.5
mplParams['mathtext.default'] = "regular"


def level1_conc_plot(level1, conc_type="Number", ucass_number=1, y_axis="Altitude"):

    if conc_type == "Number":
        unit = r'$cm^{-3}$'
        if "CYISUAData" in str(type(level1)):
            if ucass_number == 1:
                conc = level1.number_concentration1 / 1e6
            elif ucass_number == 2:
                conc = level1.number_concentration2 / 1e6
            else:
                raise ValueError("ERROR: Invalid UCASS number, only 1 or 2 accepted")
        else:
            conc = level1.number_concentration / 1e6

    elif conc_type == "Mass":
        unit = r'$kgm^{-3}$'
        if "CYISUAData" in str(type(level1)):
            if ucass_number == 1:
                conc = level1.mass_concentration1
            elif ucass_number == 2:
                conc = level1.mass_concentration2
            else:
                raise ValueError("ERROR: Invalid UCASS number, only 1 or 2 accepted")
        else:
            conc = level1.mass_concentration
    else:
        raise ValueError("ERROR: Unrecognised conc_type, \"Number\" or \"Mass\" are accepted")
    if y_axis == "Pressure":
        alt = level1.press_hpa
    elif y_axis == "Altitude":
        alt = level1.alt
    else:
        raise ValueError("ERROR: Only \"Pressure\" or \"Altitude\" on y_axis")
    up_mask = level1.up_profile_mask
    down_mask = level1.down_profile_mask
    r, c = up_mask.shape
    alt = np.delete(alt, np.s_[r:])
    conc = np.delete(conc, np.s_[r:])

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(8.3, 15))
    ax = fig.add_axes([0.2, 0.15, 0.75, 0.7])
    title_string = "Stratified %s Concentration" % conc_type
    ax.set_title(title_string, fontsize="small")
    if y_axis == "Pressure":
        ax.set_ylabel('Pressure (hPa)', fontsize="small")
    elif y_axis == "Altitude":
        ax.set_ylabel('Altitude (m)', fontsize="small")
    ax.set_xlabel('%s Concentration (%s)' % (conc_type, unit), fontsize="small")

    window = int(common.read_setting("conc_window_size"))
    filtered_conc = np.convolve(conc, np.ones((window,)) / window, mode="same")
    up_masked_conc = filtered_conc[np.where(up_mask[:, 0] == 1)]
    up_masked_sig = conc[np.where(up_mask[:, 0] == 1)]
    down_masked_conc = filtered_conc[np.where(down_mask[:, 0] == 1)]
    down_masked_sig = conc[np.where(down_mask[:, 0] == 1)]
    up_masked_alt = alt[np.where(up_mask[:, 0] == 1)]
    down_masked_alt = alt[np.where(down_mask[:, 0] == 1)]

    patch_handles = []
    marker_style = dict(linestyle='-', marker=None, markersize=5, fillstyle='none', color='C0', linewidth=0.7)
    scatter_style = dict(marker='o', color='C0', s=1, alpha=0.5)
    ax.plot(up_masked_conc, up_masked_alt, **marker_style)
    patch_handles.append(lines.Line2D([], [], linestyle='-', marker=None, color='C0', linewidth=0.7))
    ax.scatter(up_masked_sig, up_masked_alt, **scatter_style)
    patch_handles.append(lines.Line2D([], [], marker='o', color='C0', alpha=0.5, linestyle='none'))

    marker_style["color"] = 'C1'
    scatter_style["color"] = 'C1'
    ax.plot(down_masked_conc, down_masked_alt, **marker_style)
    patch_handles.append(lines.Line2D([], [], linestyle='-', marker=None, color='C1', linewidth=0.7))
    ax.scatter(down_masked_sig, down_masked_alt, **scatter_style)
    patch_handles.append(lines.Line2D([], [], marker='o', color='C1', alpha=0.5, linestyle='none'))

    leg_labels = ["Ascent with Window={w}".format(w=window), "Ascent Raw",
                  "Descent with Window={w}".format(w=window), "Descent Raw"]
    leg = Legend(ax, patch_handles, leg_labels, frameon=False, fontsize="small")
    ax.add_artist(leg)
    if y_axis == "Pressure":
        plt.gca().invert_yaxis()
    return fig


def level1_psd_plot(level1, altitude_list, ucass_number=1):

    if "CYISUAData" in str(type(level1)):
        if ucass_number == 1:
            bins = level1.bin_centres_dp_um1
        elif ucass_number == 2:
            bins = level1.bin_centres_dp_um2
        else:
            raise ValueError("ERROR: Invalid UCASS number, only 1 or 2 accepted")
    else:
        bins = level1.bin_centres_dp_um

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(12, 8))
    ax = fig.add_axes([0.2, 0.15, 0.75, 0.7])
    title_string = "Height Resolved dn/dlog(Dp)"
    ax.set_title(title_string, fontsize="small")
    ax.set_ylabel('Normalised Concentration\n(dN/dlogDp)', fontsize="small")
    ax.set_xlabel(r'Particle Diameter ($\mu m$)', fontsize="small")

    marker_styles = ["x", "o", "+", "s", "D"]
    line_styles = ['solid', 'dotted', 'dashed', 'dashdot']
    style = dict(linestyle='-', marker='none', markersize=5, fillstyle='none', color='C0', linewidth=0.7)

    sigmas = {}
    means = {}
    for alt in altitude_list:
        rows_for_mean = level1to2.fetch_row_tolerance(altitude=alt, level1_data=level1)
        sigma, mean = level1to2.mean_dn_dlogdp(level1, rows_for_mean, ucass_number=ucass_number)
        sigmas[alt] = sigma
        means[alt] = mean

    line_handles = []
    patch_handles = []
    leg_labels = []
    index = 0
    for key in means:

        while True:
            try:
                style['color'] = "C" + str(index)
                style['marker'] = marker_styles[index]
                style['linestyle'] = line_styles[index]
                break
            except IndexError:
                index -= len(line_styles)

        line_handle = ax.errorbar(bins, means[key], yerr=sigmas[key], **style)
        line_handles.append(line_handle)

        patch_handle = lines.Line2D([], [], **style)
        patch_handles.append(patch_handle)

        leg_label = key
        leg_labels.append(leg_label)

        index += 1

    leg = Legend(ax, patch_handles, leg_labels, frameon=False, fontsize="small")
    ax.add_artist(leg)
    return fig


def level1_stratified_size(level1, thresholds_list, ucass_number=1):

    if "CYISUAData" in str(type(level1)):
        if ucass_number == 1:
            bins = level1.bin_bounds_dp_um1
            dn_dlogdp = level1.dn_dlogdp1
        elif ucass_number == 2:
            bins = level1.bin_bounds_dp_um2
            dn_dlogdp = level1.dn_dlogdp2
        else:
            raise ValueError("ERROR: Invalid UCASS number, only 1 or 2 accepted")
    else:
        bins = level1.bin_bounds_dp_um
        dn_dlogdp = level1.dn_dlogdp

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(8.3, 15))
    ax = fig.add_axes([0.2, 0.15, 0.75, 0.7])
    title_string = "Particle Number Concentration - Size Resolved"
    ax.set_title(title_string, fontsize="small")
    ax.set_ylabel('Normalised Concentration\n(dN/dlogDp)', fontsize="small")
    ax.set_xlabel(r'Particle Diameter ($\mu m$)', fontsize="small")

    marker_styles = ["x", "o", "+", "s", "D"]
    line_styles = ['solid', 'dotted', 'dashed', 'dashdot']
    style = dict(linestyle='-', marker='none', markersize=5, fillstyle='none', color='C0', linewidth=0.7)

    alt, dn = level1to2.extract_dn_columns(dn_dlogdp, bins, thresholds_list)

    handles = []
    patch_handles = []
    leg_labels = []
    for i in range(dn.shape[1]):
        style['color'] = "C" + str(i)
        style['marker'] = marker_styles[i]
        style['linestyle'] = line_styles[i]

        handle = ax.plot(alt, dn, **style)
        handles.append(handle)

        patch_handle = lines.Line2D([], [], **style)
        patch_handles.append(patch_handle)

        leg_label = "%sum-%sum" % (str(thresholds_list[i]), str(thresholds_list[i+1]))
        leg_labels.append(leg_label)

    leg = Legend(ax, patch_handles, leg_labels, frameon=False, fontsize="small")
    ax.add_artist(leg)
    return fig
