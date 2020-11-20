from AirborneParticleAnalysis import level1to2
from AirborneParticleAnalysis import common
from matplotlib import pyplot as plt
from matplotlib import font_manager as plt_fnt
from matplotlib import rcParams as mplParams
import numpy as np


def _make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def _dy_dx(y_arr, x_arr):
    out = np.zeros(np.shape(y_arr))
    for x1, x2, y1, y2, i in \
            zip(x_arr[0:-1], x_arr[1:], y_arr[0:-1], y_arr[1:], [i+1 for i in (range(y_arr.shape[0]-1))]):
        out[i] = (y2 - y1) / (x2 - x1)
    return out


if __name__ == "__main__":

    prop = plt_fnt.FontProperties(family=['serif'])
    mplParams["font.family"] = prop.get_name()
    mplParams['hatch.linewidth'] = 0.5
    mplParams['mathtext.default'] = "regular"

    sam_data = level1to2.import_level1("C:\\Users\\jg17acv\\University of Hertfordshire\\AeroSAM - Documents\\Data"
                                       "\\2020\\20-09-28\\level_1\\FMITalon_PID-RM-003_20200928_10420556_00.pdat")

    mask = sam_data.up_profile_mask
    # mask = 1

    val, v_val, aoa = level1to2.check_valid_fixedwing(sam_data, 200, aoa_lim_deg=10)
    m_tof_us = np.true_divide(sam_data.m_tof, 3.0)[:, 0]
    v_m_tof_ms = np.multiply(np.true_divide(80, m_tof_us), 4.7)
    v = sam_data.vz_cms
    t = sam_data.time[:v.shape[0]]
    h = sam_data.press_hpa
    p = sam_data.pitch
    y = sam_data.yaw
    gs = sam_data.v_gnd_cms
    conc = sam_data.number_concentration / 1e6
    ws_v_arr = -1 * _dy_dx(sam_data.alt, sam_data.time)
    for prof, i in zip(mask, range(mask.shape[0])):
        if prof == 0:
            v_m_tof_ms[i] = np.nan
            conc[i] = np.nan
            val[i] = np.nan
            v_val[i] = np.nan
            aoa[i] = np.nan
            v[i] = np.nan
            t[i] = np.nan
            h[i] = np.nan
            p[i] = np.nan
            y[i] = np.nan
            gs[i] = np.nan
            ws_v_arr[i] = np.nan

    val = val[~np.isnan(val)]
    v_val = v_val[~np.isnan(v_val)]
    aoa = aoa[~np.isnan(aoa)]
    v = v[~np.isnan(v)]
    t = t[~np.isnan(t)]
    h = h[~np.isnan(h)]
    p = p[~np.isnan(p)]
    y = y[~np.isnan(y)]
    gs = gs[~np.isnan(gs)]
    ws_v_arr = ws_v_arr[~np.isnan(ws_v_arr)]
    v_m_tof_ms = v_m_tof_ms[~np.isnan(v_m_tof_ms)]
    conc = conc[~np.isnan(conc)]

    fig, ax1 = plt.subplots(1, 1, figsize=common.cm_to_inch(12, 24))
    ax1.set_ylabel('Pressure (hPa)')

    color = 'tab:red'
    ax1.set_xlabel('Number Concentration (cm-1)', color=color)
    h1 = ax1.plot(conc, h, color=color, linewidth=0.3)

    color = 'tab:blue'
    ax2 = ax1.twiny()
    h2 = ax2.plot(v / 100.0, h, color=color, linewidth=0.3)
    ax2.set_xlabel('Airspeed (m/s)', color=color)

    color = 'tab:green'
    ax3 = ax1.twiny()
    ax3.set_xlim(ax2.get_xlim())
    ax3.spines["top"].set_position(("axes", 1.1))
    _make_patch_spines_invisible(ax3)
    ax3.spines["top"].set_visible(True)
    h3 = ax3.plot(level1to2.adjust_airspeed_mtof(sam_data, "up", window=5), h, color=color, linewidth=0.3)
    ax3.set_xlabel('Adjusted Airspeed (m/s)', color=color)

    [x1, x2] = ax1.get_xlim()
    y1 = h[0]
    y2 = None
    state = val[0] * v_val[0]
    for valid, press, v_valid in zip(val, h, v_val):
        if valid == state:
            continue
        else:
            y2 = press
            if state == 1:
                color = 'tab:green'
            else:
                color = 'tab:red'
            ax1.fill_between([x1, x2], [y1, y1], [y2, y2], color=color, alpha=0.2, linewidth=0.0)
            state = valid * v_valid
            y1 = press
    y2 = h[-1]
    if state == 1:
        color = 'tab:green'
    else:
        color = 'tab:red'
    ax1.fill_between([x1, x2], [y1, y1], [y2, y2], color=color, alpha=0.2, linewidth=0.0)

    cloud_index = level1to2.detect_cloud_limits(sam_data, "up")
    ax1.plot([x1, x2], [h[cloud_index[0]], h[cloud_index[0]]])
    ax1.plot([x1, x2], [h[cloud_index[1]], h[cloud_index[1]]])

    fig.tight_layout()
    plt.gca().invert_yaxis()
    plt.show()
