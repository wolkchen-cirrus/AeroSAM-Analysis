from matplotlib import pyplot as plt
from matplotlib import font_manager as plt_fnt
from matplotlib import rcParams as mplParams
from matplotlib import lines
from matplotlib.legend import Legend
import matplotlib.ticker as plticker
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np


plt.style.use("ggplot")
prop = plt_fnt.FontProperties(family=['serif'])
mplParams["font.family"] = prop.get_name()
mplParams['hatch.linewidth'] = 0.5
mplParams['mathtext.default'] = "regular"


def _cm_to_inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(k/inch for k in tupl[0])
    else:
        return tuple(k/inch for k in tupl)


def plot_pace_dn_dlogdp(sam_array, cas_array, fssp_array=None):

    fig = plt.figure()
    fig.set_size_inches(_cm_to_inch(12, 8))
    ax = fig.add_axes([0.15, 0.2, 0.8, 0.7])
    title_string = "PSD for the SUA Mounted Instrument and Reference Instrumentation"
    ax.set_title(title_string, fontsize="small")

    line_style = ['solid', 'dotted', 'dashed', 'dashdot']
    marker_style = dict(linestyle=':', marker='x', markersize=5, fillstyle='none', color=(0, 0, 0))
    legend1_style = dict(marker='x', color=(0, 0, 0), linestyle='None', fillstyle='none')

    data_dict = {"UH-AeroSAM": sam_array, "CAS": cas_array}
    if fssp_array is not None:
        data_dict["FSSP"] = fssp_array

    patch1_handles = []
    line_handles = []
    index = 0

    for data in data_dict:

        legend1_style['linestyle'] = line_style[index]
        patch1_handle = lines.Line2D([], [], **legend1_style)
        patch1_handles.append(patch1_handle)

        marker_style['linestyle'] = line_style[index]
        line_handle = ax.plot(data[0, :], data[1, :], **marker_style)
        line_handles.append(line_handle)

        index += 1

    ax.set_ylabel("Mass Flux Ratio\nThrough Radial Region", fontsize="small")
    ax.set_xlabel("Region Elevation Above Propellers (m)", fontsize="small")
    ax.set_ylim(ymin=0)

    plt.show()
    return
