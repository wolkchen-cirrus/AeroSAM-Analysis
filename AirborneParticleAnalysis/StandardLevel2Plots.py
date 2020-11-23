from matplotlib import pyplot as plt
from matplotlib import font_manager as plt_fnt
from matplotlib import rcParams as mplParams
from matplotlib import lines
from matplotlib.legend import Legend
import numpy as np
from AirborneParticleAnalysis import common
from AirborneParticleAnalysis import StandardLevel1Plots
import scipy.stats as st
import matplotlib.gridspec as gridspec


prop = plt_fnt.FontProperties(family=['serif'])
mplParams["font.family"] = prop.get_name()
mplParams['hatch.linewidth'] = 0.5
mplParams['mathtext.default'] = "regular"
plt.rcParams['xtick.labelsize'] = "small"
plt.rcParams['ytick.labelsize'] = "small"


def level2_conc_plot(level2, prof_num=1, conc_type="Number", ucass_number=1, y_axis="Altitude", dn=None,
                     asp_lim=None, conc_lim=None, dn_lim=None, alt_lim=None):

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

    fig = plt.figure(constrained_layout=False)
    ax3 = None
    if dn is None:
        gs1 = gridspec.GridSpec(nrows=4, ncols=4, figure=fig, wspace=0)
        ax1 = fig.add_subplot(gs1[:, :2])
        ax2 = fig.add_subplot(gs1[:, 2:])
    else:
        gs = gridspec.GridSpec(ncols=4, nrows=4, figure=fig, left=0.165, right=0.98)
        gs1 = gridspec.GridSpec(nrows=4, ncols=4, figure=fig, wspace=0)
        ax1 = fig.add_subplot(gs1[:, 0])
        ax2 = fig.add_subplot(gs1[:, 1])
        plt.setp(ax2.get_yticklabels(), visible=False)
        ax3 = fig.add_subplot(gs[0:2, 2:])
        StandardLevel1Plots.level1_psd_plot(level2, dn, axis=ax3, prof_num=prof_num)

    ax4 = ax1.twiny()
    ax5 = ax2.twiny()
    ax4.set_xlabel("$Airspeed (ms^{-1})$", fontsize="small")
    ax5.set_xlabel("$Airspeed (ms^{-1})$", fontsize="small")
    fig.set_size_inches(common.cm_to_inch(18, 15))
    dt = level2.datetime
    title_datetime = dt[0:4] + "/" + dt[4:6] + "/" + dt[6:8] + " " + dt[8:10] + ":" + dt[10:12] + ":" + dt[12:14]
    title_string = "%s Stratified %s Concentration" % (title_datetime, conc_type)
    fig.suptitle(title_string)
    if y_axis == "Pressure":
        ax1.set_ylabel('Pressure (hPa)', fontsize="small")
    elif y_axis == "Altitude":
        ax1.set_ylabel('Altitude (m)', fontsize="small")
    ax1.set_xlabel('Ascending %s\nConcentration (%s)' % (conc_type, unit), fontsize="small")
    ax2.set_xlabel('Descending %s\nConcentration (%s)' % (conc_type, unit), fontsize="small")

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
    marker_style = dict(linestyle='-', marker=None, markersize=5, fillstyle='none', color='tab:red', linewidth=0.7)
    scatter_style = dict(marker='o', color='tab:red', s=1, alpha=0.5)
    ax1.plot(up_masked_conc, up_masked_alt, **marker_style)
    patch_handles.append(lines.Line2D([], [], linestyle='-', marker=None, color='tab:red', linewidth=0.7))
    ax1.scatter(up_masked_sig, up_masked_alt, **scatter_style)
    patch_handles.append(lines.Line2D([], [], marker='o', color='tab:red', alpha=0.5, linestyle='none'))

    ax2.plot(down_masked_conc, down_masked_alt, **marker_style)
    ax2.scatter(down_masked_sig, down_masked_alt, **scatter_style)

    ax4.plot(up_masked_arsp / 100.0, up_masked_alt, color=[0, 0, 0], linestyle='-', linewidth=0.3)
    ax5.plot(down_masked_arsp / 100.0, down_masked_alt, color=[0, 0, 0], linewidth=0.3, linestyle='-')

    ax1.set_xlim(xmin=0)
    ax2.set_xlim(xmin=0)
    ax4.set_xlim(xmin=0)
    ax5.set_xlim(xmin=0)
    if asp_lim is None:
        pass
    else:
        ax4.set_xlim(xmin=asp_lim[0], xmax=asp_lim[1])
        ax5.set_xlim(xmin=asp_lim[0], xmax=asp_lim[1])
    if alt_lim is None:
        pass
    else:
        ax1.set_ylim(ymin=alt_lim[0], ymax=alt_lim[1])
        ax2.set_ylim(ymin=alt_lim[0], ymax=alt_lim[1])
    if conc_lim is None:
        pass
    else:
        ax1.set_xlim(xmin=conc_lim[0], xmax=conc_lim[1])
        ax2.set_xlim(xmin=conc_lim[0], xmax=conc_lim[1])
    if dn_lim is None:
        pass
    elif dn is not None:
        ax3.set_ylim(ymin=dn_lim[0], ymax=dn_lim[1])
    else:
        pass

    aoa_mask = level2.aoa_mask
    vz_mask = level2.vz_mask

    color = [0, 0, 0]
    profiles = [up_mask, down_mask]
    types = [vz_mask, aoa_mask]
    axes = [ax4, ax5]
    alphas = [0.5, 0.2]
    ax6 = fig.add_subplot(gs1[-1, -1])
    dh_1 = ax6.fill_between([], [], [], color=color, alpha=alphas[0], linewidth=0.3, linestyle='-')
    dh_2 = ax6.fill_between([], [], [], color=color, alpha=alphas[1], linewidth=0.3, linestyle='-')
    leg_handles = [dh_1, dh_2, patch_handles[0], patch_handles[1]]
    leg_names = ["Airspeed > 20", "Angle of Attack > 10",
                 "Filtered Concentration\nwith Window Size of 15", "Raw Concentration"]
    leg = Legend(ax6, leg_handles, leg_names, frameon=False, fontsize="small")
    ax6.add_artist(leg)
    ax6.axis('off')
    for mask, axis in zip(profiles, axes):
        masked_alt = alt[np.where(mask[:, prof_num-1] == 1)]
        [x1, x2] = axis.get_xlim()
        for val_type, a in zip(types, alphas):
            y1 = masked_alt[0]
            masked_val_type = val_type[np.where(mask[:, prof_num - 1] == 1)]
            q_state = masked_val_type[0]
            for q_valid, press in zip(masked_val_type, masked_alt):
                if q_valid == q_state:
                    continue
                else:
                    y2 = press
                    if q_state == 1:
                        alpha = 0
                    else:
                        alpha = a
                    axis.fill_between([x1, x2], [y1, y1], [y2, y2], color=color, alpha=alpha,
                                      linewidth=0.3, linestyle='-')
                    q_state = q_valid
                    y1 = press

    if y_axis == "Pressure":
        plt.gca().invert_yaxis()

    tick1_remove_list = [ax2, ax5]
    for ax in tick1_remove_list:
        x_ticks = ax.xaxis.get_major_ticks()
        x_ticks[0].label1.set_visible(False)
        x_ticks[0].label.set_visible(False)
        x_ticks[0].label2.set_visible(False)

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
