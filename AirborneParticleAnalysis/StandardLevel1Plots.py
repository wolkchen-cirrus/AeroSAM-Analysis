from matplotlib import pyplot as plt
from matplotlib import font_manager as plt_fnt
from matplotlib import rcParams as mplParams
from matplotlib import lines
from matplotlib.legend import Legend
# import matplotlib.ticker as plticker
from AirborneParticleAnalysis import common
from AirborneParticleAnalysis import level1to2


# plt.style.use("ggplot")
prop = plt_fnt.FontProperties(family=['serif'])
mplParams["font.family"] = prop.get_name()
mplParams['hatch.linewidth'] = 0.5
mplParams['mathtext.default'] = "regular"


def level1_conc_plot(conc, alt, conc_sigma, conc_type="Number"):

    if conc_type == "Number":
        unit = r'$cm^{-1}$'
    elif conc_type == "Mass":
        unit = r'$gcm^{-1}$'
    else:
        raise ValueError("ERROR: Unrecognised conc_type, \"Number\" or \"Mass\" are accepted")

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(8.3, 10))
    ax = fig.add_axes([0.2, 0.15, 0.75, 0.7])
    title_string = "Stratified %s Concentration" % conc_type
    ax.set_title(title_string, fontsize="small")
    ax.set_ylabel('Altitude (m)', fontsize="small")
    ax.set_xlabel('%s Concentration (%s)' % conc_type, unit, fontsize="small")

    marker_style = dict(linestyle='-', marker='none', markersize=5, fillstyle='none', color=(0, 0, 0), linewidth=0.7)
    ax.errorbar(conc, alt, yerr=conc_sigma, **marker_style)

    return


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

        patch_handle = lines.Line2D([], [], **style)
        patch_handles.append(patch_handle)

        style['color'] = "C" + str(index)
        style['marker'] = marker_styles[index]
        style['linestyle'] = line_styles[index]
        line_handle = ax.errorbar(bins, means[key], yerr=sigmas[key], **style)
        line_handles.append(line_handle)

        leg_label = key
        leg_labels.append(leg_label)

        index += 1

    leg = Legend(ax, patch_handles, leg_labels, frameon=False, fontsize="small")
    ax.add_artist(leg)

    return


def level1_stratified_size(level1, bin_list, ucass_number=1):

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
    fig.set_size_inches(common.cm_to_inch(12, 8))
    ax = fig.add_axes([0.2, 0.15, 0.75, 0.7])
    title_string = "Particle Number Concentration - Size Resolved"
    ax.set_title(title_string, fontsize="small")
    ax.set_ylabel('Normalised Concentration\n(dN/dlogDp)', fontsize="small")
    ax.set_xlabel(r'Particle Diameter ($\mu m$)', fontsize="small")

    marker_styles = ["x", "o", "+", "s", "D"]
    line_styles = ['solid', 'dotted', 'dashed', 'dashdot']
    style = dict(linestyle='-', marker='none', markersize=5, fillstyle='none', color='C0', linewidth=0.7)

    return
