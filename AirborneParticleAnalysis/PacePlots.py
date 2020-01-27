from matplotlib import pyplot as plt
from matplotlib import font_manager as plt_fnt
from matplotlib import rcParams as mplParams
from matplotlib import lines
from matplotlib.legend import Legend
import matplotlib.ticker as plticker
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
from AirborneParticleAnalysis import common


plt.style.use("ggplot")
prop = plt_fnt.FontProperties(family=['serif'])
mplParams["font.family"] = prop.get_name()
mplParams['hatch.linewidth'] = 0.5
mplParams['mathtext.default'] = "regular"


def plot_rebin_1to1(data_ref, data_sam, lr_x, lr_y):

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(12, 8))
    ax = fig.add_axes([0.15, 0.2, 0.8, 0.7])
    title_string = "PSD for the SUA Mounted Instrument and Reference Instrumentation"
    ax.set_title(title_string, fontsize="small")

    marker_list = ['x', '+', 'x', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']
    marker_style = dict(linestyle='none', marker='x', markersize=5, fillstyle='none', color=(0, 0, 0))
    legend1_style = dict(marker='x', color=(0, 0, 0), linestyle='None', fillstyle='none')

    return


def plot_pace_dn_dlogdp(data_dict, sam_bins=None, cas_bins=None, fssp_bins=None):

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(12, 8))
    ax = fig.add_axes([0.15, 0.2, 0.8, 0.7])
    title_string = "PSD for the SUA Mounted Instrument and Reference Instrumentation"
    ax.set_title(title_string, fontsize="small")

    line_style = ['solid', 'dotted', 'dashed', 'dashdot']
    marker_style = dict(linestyle=':', marker='x', markersize=5, fillstyle='none', color=(0, 0, 0))
    legend1_style = dict(marker='x', color=(0, 0, 0), linestyle='None', fillstyle='none')

    patch1_handles = []
    line_handles = []
    leg_labels = []
    index = 0
    bins = None
    for key in data_dict:

        data = data_dict[key]

        if "SAM" in key:
            bins = sam_bins
            legend1_style['linestyle'] = line_style[index]
            patch1_handle = lines.Line2D([], [], **legend1_style)
            patch1_handles.append(patch1_handle)

            marker_style['linestyle'] = line_style[index]
            line_handle = ax.plot(bins, data, **marker_style)
            line_handles.append(line_handle)

            leg_label = "SUA with UCASS-V2"
            leg_labels.append(leg_label)

            index += 1

    if bins is None:
        raise ValueError("ERROR: No bins specified")

    bins = None
    for key in data_dict:

        data = data_dict[key]
        leg_label = None

        if ("CAS" in key) or ("FSSP" in key):
            if "CAS" in key:
                bins = cas_bins
                leg_label = "CAS at %sm" % str(int(common.read_setting("station_altitude_asl_mm"))/1000)
            elif "FSSP" in key:
                bins = fssp_bins
                leg_label = "FSSP at %sm" % common.read_setting("station_altitude_asl_mm")
            legend1_style['linestyle'] = line_style[index]
            patch1_handle = lines.Line2D([], [], **legend1_style)
            patch1_handles.append(patch1_handle)

            marker_style['linestyle'] = line_style[index]
            line_handle = ax.plot(bins, data, **marker_style)
            line_handles.append(line_handle)

            leg_labels.append(leg_label)

            index += 1

    if bins is None:
        raise ValueError("ERROR: No bins specified")

    ax.set_ylabel('Normalised Concentration\n(dN/dlogDp)', fontsize="small")
    ax.set_xlabel(r'Particle Diameter ($\mu m$)', fontsize="small")
    ax.set_ylim(ymin=0)
    ax.set_xlim(xmin=0)

    leg = Legend(ax, patch1_handles, leg_labels)
    ax.add_artist(leg)

    plt.show()
    return
