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


def level2_conc_plot(level2, prof_num=1, conc_type="Number", ucass_number=1, y_axis="Altitude"):

    if conc_type == "Number":
        unit = r'$cm^{-3}$'
        if "CYISUAData" in str(type(level2)):
            if ucass_number == 1:
                conc = level2.number_concentration1 / 1e6
            elif ucass_number == 2:
                conc = level2.number_concentration2 / 1e6
            else:
                raise ValueError("ERROR: Invalid UCASS number, only 1 or 2 accepted")
        else:
            conc = level2.number_concentration / 1e6

    elif conc_type == "Mass":
        unit = r'$kgm^{-3}$'
        if "CYISUAData" in str(type(level2)):
            if ucass_number == 1:
                conc = level2.mass_concentration1
            elif ucass_number == 2:
                conc = level2.mass_concentration2
            else:
                raise ValueError("ERROR: Invalid UCASS number, only 1 or 2 accepted")
        else:
            conc = level2.mass_concentration
    else:
        raise ValueError("ERROR: Unrecognised conc_type, \"Number\" or \"Mass\" are accepted")
    if y_axis == "Pressure":
        alt = level2.press_hpa
    elif y_axis == "Altitude":
        alt = level2.alt
    else:
        raise ValueError("ERROR: Only \"Pressure\" or \"Altitude\" on y_axis")
    up_mask = level2.up_profile_mask
    down_mask = level2.down_profile_mask
    arsp = level2.adjusted_airspeed
    r, c = up_mask.shape
    alt = np.delete(alt, np.s_[r:])
    conc = np.delete(conc, np.s_[r:])

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(8.3, 15))
    ax = fig.add_axes([0.2, 0.15, 0.75, 0.7])
    ax2 = ax.twiny()
    ax2.set_xlabel("$Airspeed (ms^{-1})$")
    dt = level2.datetime
    title_datetime = dt[0:4] + "/" + dt[4:6] + "/" + dt[6:8] + " " + dt[8:10] + ":" + dt[10:12] + ":" + dt[12:14]
    title_string = "%s\nStratified %s Concentration" % (title_datetime, conc_type)
    ax.set_title(title_string, fontsize="small")
    if y_axis == "Pressure":
        ax.set_ylabel('Pressure (hPa)', fontsize="small")
    elif y_axis == "Altitude":
        ax.set_ylabel('Altitude (m)', fontsize="small")
    ax.set_xlabel('%s Concentration (%s)' % (conc_type, unit), fontsize="small")

    window = int(common.read_setting("conc_window_size"))
    filtered_conc = np.convolve(conc, np.ones((window,)) / window, mode="same")
    up_masked_conc = filtered_conc[np.where(up_mask[:, prof_num-1] == 1)]
    up_masked_sig = conc[np.where(up_mask[:, prof_num-1] == 1)]
    down_masked_conc = filtered_conc[np.where(down_mask[:, prof_num-1] == 1)]
    down_masked_sig = conc[np.where(down_mask[:, prof_num-1] == 1)]
    up_masked_alt = alt[np.where(up_mask[:, prof_num-1] == 1)]
    down_masked_alt = alt[np.where(down_mask[:, prof_num-1] == 1)]
    up_masked_arsp = arsp[np.where(up_mask[:, prof_num-1] == 1)]
    down_masked_arsp = arsp[np.where(down_mask[:, prof_num - 1] == 1)]

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

    ax2.plot(up_masked_arsp / 100.0, up_masked_alt, color=[0, 0, 0], linestyle='-', linewidth=0.3)
    ax2.plot(down_masked_arsp / 100.0, down_masked_alt, color=[0, 0, 0], linewidth=0.3, linestyle='-.')

    leg_labels = ["Ascent with Window={w}".format(w=window), "Ascent Raw",
                  "Descent with Window={w}".format(w=window), "Descent Raw"]
    leg = Legend(ax, patch_handles, leg_labels, frameon=False, fontsize="small")
    ax.add_artist(leg)

    ax.set_xlim(xmin=0)

    aoa_mask = level2.aoa_mask
    vz_mask = level2.vz_mask

    [x1, x2] = ax.get_xlim()
    color = [0, 0, 0]
    y1 = alt[0]
    v_state = vz_mask[0]
    a_state = aoa_mask[0]
    for a_valid, press in zip(aoa_mask, alt):
        if a_valid == a_state:
            continue
        else:
            y2 = press
            if a_state == 1:
                alpha = 0
            else:
                alpha = 0.2
            ax.fill_between([x1, x2], [y1, y1], [y2, y2], color=color, alpha=alpha, linewidth=0.0)
            a_state = a_valid
            y1 = press
    for v_valid, press in zip(vz_mask, alt):
        if v_valid == v_state:
            continue
        else:
            y2 = press
            if v_state == 1:
                alpha = 0
            else:
                alpha = 0.5
            ax.fill_between([x1, x2], [y1, y1], [y2, y2], color=color, alpha=alpha, linewidth=0.0)
            v_state = v_valid
            y1 = press

    if y_axis == "Pressure":
        plt.gca().invert_yaxis()

    return fig, title_string


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

