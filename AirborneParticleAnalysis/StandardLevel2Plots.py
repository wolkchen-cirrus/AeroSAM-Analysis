from matplotlib import pyplot as plt
from matplotlib import font_manager as plt_fnt
from matplotlib import rcParams as mplParams
from matplotlib import lines
from matplotlib.legend import Legend
import numpy as np
from AirborneParticleAnalysis import common
import scipy.stats as st


prop = plt_fnt.FontProperties(family=['serif'])
mplParams["font.family"] = prop.get_name()
mplParams['hatch.linewidth'] = 0.5
mplParams['mathtext.default'] = "regular"


def plot_rebin_1to1(data_1, data_2, bin_centres, regression_data, names, contour=None):

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(12, 12))
    ax = fig.add_axes([0.2, 0.15, 0.75, 0.7])
    title_string = names[2] + " " + "Integrated dn/dlog(Dp)"
    ax.set_title(title_string, fontsize="small")
    ax.set_ylabel(names[1], fontsize="small")
    ax.set_xlabel(names[0], fontsize="small")

    x12 = regression_data[0]
    y12 = regression_data[1]
    r2 = regression_data[2]
    m = regression_data[5]

    data_z = []
    index = 0
    while True:
        if len(data_z) == len(data_1):
            break
        data_z.append(bin_centres[index])
        index += 1
        if index == 16:
            index = 0

    bar = plt.cm.get_cmap('RdYlBu')
    bar.label = r'Bin Centre \mu m'

    marker_style = dict(linestyle='none', marker='x', markersize=5, fillstyle='none', linewidth=0.7)
    handle = ax.scatter(data_1, data_2, c=data_z, cmap=bar, s=0.3)
    plt.colorbar(handle)

    if contour is None:
        pass
    else:
        delta_x = (max(data_1) - min(data_1)) / 10
        delta_y = (max(data_2) - min(data_2)) / 10
        x_min = min(data_1) - delta_x
        x_max = max(data_1) + delta_x
        y_min = min(data_2) - delta_y
        y_max = max(data_2) + delta_y
        xx, yy = np.mgrid[x_min:x_max:100j, y_min:y_max:100j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([data_1, data_2])
        kernel = st.gaussian_kde(values)
        f = np.reshape(kernel(positions).T, xx.shape)
        ax.imshow(np.rot90(f), cmap='gray_r', extent=[x_min, x_max, y_min, y_max])

    marker_style["color"] = (0, 0, 0)

    marker_style["linestyle"] = '-.'
    marker_style["marker"] = ''
    ax.plot(x12, y12, **marker_style)
    h1 = lines.Line2D([], [], **marker_style)

    marker_style["linestyle"] = 'solid'
    ax.plot([0, x12[-1]], [0, x12[-1]], **marker_style)
    h2 = lines.Line2D([], [], **marker_style)

    ax.text(0.2, 0.8, "$r^{2} = $%f\n$m = $%f" % (r2, m),
            fontsize="small", ha='center', va='center', transform=ax.transAxes)

    ax.set_ylim(ymin=0)
    ax.set_xlim(xmin=0)

    leg = Legend(ax, [h1, h2], ["Regression Line", "y = x"], frameon=False, fontsize="small", loc=2)
    ax.add_artist(leg)

    return fig, title_string


def plot_cint_dn(level1, end_swipe, resolve_type="mean"):

    if "CYISUAData" in str(type(level1)):
        bins = [level1.bin_centres_dp_um1, level1.bin_centres_dp_um2]
        dn_dlogdp = [level1.dn_dlogdp1.values(), level1.dn_dlogdp2.values()]
    else:
        bins = [level1.bin_centres_dp_um]
        dn_dlogdp = [level1.dn_dlogdp]

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(12, 8))
    ax = fig.add_axes([0.2, 0.15, 0.75, 0.7])
    dt = level1.datetime
    title_datetime = dt[0:4] + "/" + dt[4:6] + "/" + dt[6:8] + " " + dt[8:10] + ":" + dt[10:12] + ":" + dt[12:14]
    title_string = "%s\nColumn-Mean dn/dlog(Dp)" % title_datetime
    ax.set_title(title_string, fontsize="small")
    ax.set_ylabel('dN/dlogDp', fontsize="small")
    ax.set_xlabel(r'Particle Diameter ($\mu m$)', fontsize="small")

    dn_dlogdp_arr = []
    for dn in dn_dlogdp:
        del dn[0:end_swipe]
        del dn[-end_swipe:-1]
        dn_arr = np.array([np.array(dni) for dni in dn])
        dn_dlogdp_arr.append(dn_arr)

    dn_p = []
    for dn in dn_dlogdp:
        if resolve_type == "mean":
            dn_p.append(np.mean(dn, axis=0))
        if resolve_type == "integrated":
            dn_p.append(np.sum(dn, axis=0))

    marker_styles = ["x", "o", "+", "s", "D"]
    line_styles = ['solid', 'dotted', 'dashed', 'dashdot']
    style = dict(linestyle='-', marker='none', markersize=5, fillstyle='none', color='C0', linewidth=0.7)

    handles = []
    patch_handles = []
    leg_labels = []
    i = 0
    for dn in dn_p:
        style['color'] = "C" + str(i)
        style['marker'] = marker_styles[i]
        style['linestyle'] = line_styles[i]

        handle = ax.plot(bins[i], dn, **style)
        handles.append(handle)

        patch_handle = lines.Line2D([], [], **style)
        patch_handles.append(patch_handle)

        leg_label = "UCASS %s" % str(i+1)
        leg_labels.append(leg_label)

        i += 1

    leg = Legend(ax, patch_handles, leg_labels, frameon=False, fontsize="small")
    ax.add_artist(leg)

    ax.set_ylim(ymin=0)
    ax.set_xlim(xmin=0)

    return fig, title_string

