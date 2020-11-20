import numpy as np
import cPickle as Pickle
from AirborneParticleAnalysis import common
from scipy import stats
from scipy.signal import find_peaks
from scipy.interpolate import interp1d


def import_level1(level1_path):
    with open(level1_path) as level1_file:
        level1_data = Pickle.load(level1_file)
    return level1_data


def fetch_row(altitude=None, time=None, level1_data=None, profile="Up"):

    if "SUAData" in str(type(level1_data)):
        try:
            key_col = level1_data.alt
            if not isinstance(key_col, np.ndarray):
                try:
                    key_col = np.asarray(key_col)
                except(ValueError, TypeError):
                    raise TypeError("ERROR: Incompatible data type")
            row_value = altitude
            if profile == "Up":
                prof_mask = level1_data.up_profile_mask
            elif profile == "Down":
                prof_mask = level1_data.down_profile_mask
            else:
                raise ValueError("ERROR: Profile is either \"Up\" or \"Down\" as str")
            (r, cols) = prof_mask.shape
            rows = []
            if "CYI" in str(type(level1_data)):
                for i in range(cols):
                    key_col = key_col[np.where(prof_mask[:, 0] == 1)]
                    diff_col = abs(key_col - row_value)
                    min_diff = np.amin(diff_col)
                    min_diff_index = list(np.where(diff_col == min_diff)[0])[0]
                    buf = key_col[min_diff_index][0]
                    rows.append(buf)
            else:
                for i in range(cols):
                    diff_col = np.multiply(abs(key_col - row_value) - 100000000, np.reshape(prof_mask[:, i], (r, 1)))
                    min_diff = np.amin(diff_col)
                    min_diff_index = np.where(diff_col == min_diff)
                    buf = key_col[min_diff_index[0][0]][0]
                    rows.append(buf)
        except AttributeError:
            raise AttributeError("ERROR: level1_data object problem")

    elif ("StaticCASData" in str(type(level1_data))) or ("StaticFSSPData" in str(type(level1_data))):
        try:
            key_col = level1_data.time
            if not isinstance(key_col, np.ndarray):
                try:
                    key_col = np.asarray(key_col)
                except(ValueError, TypeError):
                    raise TypeError("ERROR: Incompatible data type")
            row_value = time
            rows = []
            diff_col = abs(key_col - row_value)
            min_diff = np.amin(diff_col)
            min_diff_index = np.where(diff_col == min_diff)
            buf = key_col[min_diff_index[0][0]][0]
            rows.append(buf)
        except AttributeError:
            raise AttributeError("ERROR: level1_data object problem")

    else:
        raise ValueError("ERROR: Unrecognised data object")

    return rows


def fetch_row_tolerance(altitude=None, time=None, level1_data=None, profile="Up"):

    if "SUAData" in str(type(level1_data)):
        try:
            key_col = level1_data.alt
            if "CYI" in str(type(level1_data)):
                tol = float(common.read_setting("height_mean_tolerance_metres"))
            else:
                tol = float(common.read_setting("height_mean_tolerance_metres")) * 1000.0
            if not isinstance(key_col, np.ndarray):
                try:
                    key_col = np.asarray(key_col)
                except(ValueError, TypeError):
                    raise TypeError("ERROR: Incompatible data type")
            row_value = altitude
            if profile == "Up":
                prof_mask = level1_data.up_profile_mask
            elif profile == "Down":
                prof_mask = level1_data.down_profile_mask
            else:
                raise ValueError
            (r, cols) = prof_mask.shape
            rows = []
            if "CYI" in str(type(level1_data)):
                for i in range(cols):
                    key_col = key_col[np.where(prof_mask[:, 0] == 1)]
                    offset = np.where(np.trim_zeros(prof_mask[:, 0], 'b') == 0)[0].shape

                    diff_col_u = abs(key_col - (row_value+tol))
                    min_diff_u = np.amin(diff_col_u)
                    min_diff_index_u = list(np.where(diff_col_u == min_diff_u)[0])[0]

                    diff_col_l = abs(key_col - (row_value-tol))
                    min_diff_l = np.amin(diff_col_l)
                    min_diff_index_l = list(np.where(diff_col_l == min_diff_l)[0])[0]

                    if min_diff_index_u == min_diff_index_l:
                        buf = key_col[min_diff_index_l][0]
                    else:
                        buf = key_col[min_diff_index_l:min_diff_index_u]
                        buf = list(buf[:])
                        buf = [val[0] for val in buf]

                    rows.append((buf, range(offset + min_diff_index_l, offset + min_diff_index_u)))
            else:
                for i in range(cols):
                    diff_col_l = \
                        np.multiply(abs(key_col - (row_value-tol)) - 100000000, np.reshape(prof_mask[:, i], (r, 1)))
                    min_diff_l = np.amin(diff_col_l)
                    min_diff_index_l = np.where(diff_col_l == min_diff_l)

                    diff_col_u = \
                        np.multiply(abs(key_col - (row_value+tol)) - 100000000, np.reshape(prof_mask[:, i], (r, 1)))
                    min_diff_u = np.amin(diff_col_u)
                    min_diff_index_u = np.where(diff_col_u == min_diff_u)

                    buf = list(key_col[min_diff_index_l[0][0]:min_diff_index_u[0][0]].flatten())

                    rows.append((buf, range(min_diff_index_l[0][0], min_diff_index_u[0][0])))
        except AttributeError:
            raise AttributeError("ERROR: level1_data object problem")

    elif ("StaticCASData" in str(type(level1_data))) or ("StaticFSSPData" in str(type(level1_data))):
        try:
            key_col = level1_data.time
            tol = float(common.read_setting("time_mean_tolerance_seconds"))
            if not isinstance(key_col, np.ndarray):
                try:
                    key_col = np.asarray(key_col)
                except(ValueError, TypeError):
                    raise TypeError("ERROR: Incompatible data type")
            row_value = time
            rows = []

            diff_col_l = abs(key_col - (row_value - tol))
            min_diff_l = np.amin(diff_col_l)
            min_diff_index_l = np.where(diff_col_l == min_diff_l)

            diff_col_u = abs(key_col - (row_value + tol))
            min_diff_u = np.amin(diff_col_u)
            min_diff_index_u = np.where(diff_col_u == min_diff_u)

            buf = list(key_col[min_diff_index_l[0][0]:min_diff_index_u[0][0]].flatten())

            rows.append((buf, range(min_diff_index_l[0][0], min_diff_index_u[0][0])))
        except AttributeError:
            raise AttributeError("ERROR: level1_data object problem")

    else:
        raise ValueError("ERROR: Unrecognised data object")

    return [(i, j) for i, j in zip(rows[0][0], rows[0][1])]


def mean_dn_dlogdp(level1_data, rows, ucass_number=1):

    if isinstance(rows, list) and isinstance(rows[0], list):
        rows = rows[0]
        if isinstance(rows[0], list):
            raise ValueError("ERROR: Pass only one profile into function")

    if "CYISUAData" in str(type(level1_data)):

        if ucass_number == 1:
            dn_dlogdp = level1_data.dn_dlogdp1
        elif ucass_number == 2:
            dn_dlogdp = level1_data.dn_dlogdp2
        else:
            raise ValueError("ERROR: For CYI, only 2 UCASS' are configured")

        if not isinstance(rows, list):
            print "INFO: No list specified so returning inputs"
            return dn_dlogdp[rows]

        data_arr = np.zeros([len(rows), len(dn_dlogdp[dn_dlogdp.keys()[0]])])
        (r, c) = data_arr.shape
        for i in range(r):
            for j in range(c):
                data_arr[i, j] = dn_dlogdp[rows[i]][j]

        dn_mean = np.mean(data_arr, axis=0)
        dn_std = np.std(data_arr, axis=0)

        return [dn_mean, dn_std]

    else:

        data_arr = np.zeros([len(rows), len(level1_data.dn_dlogdp[rows[0]])])
        (r, c) = data_arr.shape
        for i in range(r):
            for j in range(c):
                data_arr[i, j] = level1_data.dn_dlogdp[rows[i]][j]

        dn_mean = np.mean(data_arr, axis=0)
        dn_std = np.std(data_arr, axis=0)

    return [dn_mean, dn_std]


def get_time_from_alt(sua_data, alt_exact):

    if ("StaticCASData" in str(type(sua_data))) or ("StaticFSSPData" in str(type(sua_data))):
        raise ValueError("ERROR: Only SAM data can be input into this function")

    alt = sua_data.alt
    row = np.where(alt == alt_exact)[0]
    start_time = common.hhmmss_to_sec(sua_data.datetime[-6:])
    time = sua_data.time

    time_offset = time[row]

    if (time_offset - time[row-1]) < 0:
        time_offset = time[row-1]

    time_offset = time_offset[0]

    row_ref = time[0]
    if (row_ref - 1000000) < 0:
        row_ref = time[1]

    time_offset_corrected = time_offset - row_ref

    return int(start_time + time_offset_corrected[0])


def rebin_dn_dlogdp(dn, bins, new_bins):

    new_dn = []
    for xq in new_bins:

        extrapolate = False
        x2 = None
        y2 = None
        dy_dx = None
        index = 0
        for x1 in bins:
            if x1 >= xq:
                break
            else:
                index += 1

        try:
            x1 = bins[index-1]
            x2 = bins[index]
            y1 = dn[index-1]
            y2 = dn[index]
        except IndexError:
            extrapolate = True
            if index == 0:
                x1 = -bins[index]
                y1 = dn[index]
                dy_dx = -(dn[index] - dn[index+1])/(bins[index] - bins[index-1])
            else:
                try:
                    x1 = bins[index-1]
                    y1 = dn[index-1]
                    dy_dx = (dn[index] - dn[index-1])/(bins[index] - bins[index-1])
                except IndexError:
                    break

        if extrapolate is True:
            yq = _linear_extrapolate(dy_dx, x1, y1, xq)
        else:
            yq = _linear_interpolate(x1, y1, x2, y2, xq)
        new_dn.append(yq)

    if len(new_dn) >= len(new_bins):
        del new_dn[-1]

    integrated_dn_steps = _discrete_integral_2d(new_bins, new_dn, bins, list(dn))

    return integrated_dn_steps


def linear_regression(x_list, y_list):
    [m, c, r, p, e] = stats.linregress(x_list, y_list)
    x12 = [0, max(x_list)]
    y12 = [c, m * max(x_list) + c]
    r2 = r ** 2
    return [x12, y12, r2, p, e, m]


def rma_regression(x_list, y_list):
    """rma_regression - http://doi.wiley.com/10.1002/9781118445112.stat07912
    :param x_list: List of x data
    :param y_list: List of y data
    """

    sx = stats.tstd(x_list)
    sy = stats.tstd(y_list)

    avg_x = stats.tmean(x_list)
    avg_y = stats.tmean(y_list)

    m = sy/sx
    c = avg_y - m * avg_x
    r = stats.pearsonr(x_list, y_list)[0]
    r2 = r ** 2

    x12 = [0, max(x_list)]
    y12 = [c, m * max(x_list) + c]

    return [x12, y12, r2, None, None, m]


def _linear_extrapolate(dy_dx, x1, y1, xq):
    yq = dy_dx * (xq - x1) + y1
    return yq


def _linear_interpolate(x1, y1, x2, y2, xq):

    dy_dx_1 = (y2-y1)/(x2-x1)
    yq = dy_dx_1*(xq-x1) + y1

    return yq


def _discrete_integral_2d(qx, qy, fx, fy):

    if (not isinstance(qx, list)) or (not isinstance(qy, list)) or \
            (not isinstance(fx, list)) or (not isinstance(fy, list)):
        raise TypeError("ERROR: All inputs into this function must be lists")

    index = 0
    total_bin_areas = []
    for x1, y1 in zip(qx, qy):

        try:
            x2 = qx[index+1]
            y2 = qy[index+1]
        except IndexError:
            break

        x_mediums = []
        y_mediums = []
        for x, y in zip(fx, fy):
            if x1 < x < x2:
                x_mediums.append(x)
                y_mediums.append(y)

        x_in_bin = [x1] + x_mediums + [x2]
        y_in_bin = [y1] + y_mediums + [y2]

        step = 0
        bin_area = 0
        for xi1, yi1 in zip(x_in_bin, y_in_bin):
            try:
                xi2 = x_in_bin[step+1]
                yi2 = y_in_bin[step+1]
            except IndexError:
                break
            hi = xi2 - xi1
            bin_area += (yi1 + yi2) / 2 * hi
            step += 1

        total_bin_areas.append(bin_area)

        index += 1

    return total_bin_areas


def extract_dn_columns(dn_dict, bins, new_bins, mask):
    new_dn_dict = {}
    index = 0
    for key in dn_dict:
        if mask[index, 0] == 0:
            index += 1
            continue
        new_dn_dict[key] = rebin_dn_dlogdp(dn_dict[key], bins, new_bins)
        index += 1
    return np.array(new_dn_dict.keys()), np.array(new_dn_dict.values())


def percent_diff(d1, d2):
    if type(d1) != type(d2):
        raise TypeError
    elif isinstance(d1, list) or isinstance(d2, list):
        pd = [100.0 * float(abs(l1 - l2) / ((l1 + l2) / 2)) for l1, l2 in zip(d1, d2)]
    else:
        pd = 100.0 * float(abs(d1 - d2) / ((d1 + d2) / 2))
    return pd


def percent_err(d1, d2):
    if type(d1) != type(d2):
        raise TypeError
    elif isinstance(d1, list) or isinstance(d2, list):
        pd = [100.0 * float(abs(l1 - l2)/l1) for l1, l2 in zip(d1, d2)]
    else:
        pd = 100.0 * float(abs(d1 - d2) / d1)
    return pd


def num_conc_range(level1, bins=None, ucass_number=1):
    sv = level1.sample_volume_m3 * 1e6
    if "CYISUAData" in str(type(level1)):
        if ucass_number == 1:
            counts = level1.raw_counts1
        elif ucass_number == 2:
            counts = level1.raw_counts2
        else:
            raise ValueError
        if bins is None:
            pass
        else:
            counts = counts[:, bins]
        size = counts.shape
        num_conc_buf = np.zeros([size[0], 1])
        for i in range(size[0]):
            if sv[i] == 0:
                num_conc_buf[i, 0] = 0
            else:
                num_conc_buf[i, 0] = float(sum(counts[i, :])) / sv[i]
    else:
        counts = level1.raw_counts
        if bins is None:
            pass
        else:
            counts = counts[:, bins]
        size = counts.shape
        num_conc_buf = np.zeros([size[0], 1])
        for i in range(size[0]):
            if sv[i] == 0:
                num_conc_buf[i, 0] = 0
            else:
                num_conc_buf[i, 0] = float(sum(counts[i, :])) / sv[i]

    return num_conc_buf


def check_valid_fixedwing(level1, wa_deg, aoa_lim_deg=10, airspeed_lim_ms=20, airspeed_type="normal"):
    try:
        _ = level1.aoa_mask
        _ = level1.vz_mask
    except(NameError, AttributeError):
        raise AttributeError("ERROR: No aoa_mask or vz_mask, data may not be FW")
    if "CYISUAData" in str(type(level1)):
        raise NotImplementedError("ERROR: Not implemented for CYI data")

    if airspeed_type == "normal":
        vz = level1.vz_cms
    elif airspeed_type == "adjusted":
        vz = level1.adjusted_vz_cms
    else:
        raise ValueError("ERROR: Invalid airspeed type")
    ws_v_arr = -1 * _dy_dx(level1.alt, level1.time)
    pitch = level1.pitch
    yaw = level1.yaw
    gs = level1.v_gnd_cms

    wa_rad = np.deg2rad(wa_deg)

    vz_mask = np.zeros(np.shape(vz))
    for val, i in zip(vz, range(vz.shape[0])):
        if val < airspeed_lim_ms * 100:
            vz_mask[i] = 1
        else:
            pass

    aoa_mask = np.zeros(np.shape(vz))
    aoa_store = np.zeros(np.shape(vz))
    for p, y, v, g, i, ws_v in zip(pitch, yaw, vz, gs, range(vz.shape[0]), ws_v_arr):
        p_rad = np.deg2rad(p[0])
        y_rad = np.deg2rad(y[0])
        v_ms = np.true_divide(v[0], 100.0)
        g_ms = np.true_divide(g[0], 100.0)
        _as = np.array([np.sin(y_rad)*np.cos(p_rad), np.cos(y_rad)*np.cos(p_rad), np.sin(p_rad)])
        _as = np.multiply(v_ms, np.divide(_as, _mag(_as)))
        _gnd_as = np.multiply(_as, np.cos(p_rad))
        _gnd_as[-1] = 0
        _gs = np.array([np.sin(y_rad), np.cos(y_rad), 0])
        _gs = np.multiply(g_ms, np.divide(_gs, _mag(_gs)))
        _ws_h = np.divide(_gnd_as - _gs, np.deg2rad(180) - np.cos(wa_rad))
        _ws = _ws_h + np.array([0, 0, ws_v[0]])
        _r = _ws + _as
        aoa = np.arccos((np.true_divide(np.dot(_r, _as), np.multiply(_mag(_r), _mag(_as)))))
        aoa_store[i] = np.rad2deg(aoa)

        if np.rad2deg(aoa) < aoa_lim_deg:
            aoa_mask[i] = 1
        else:
            pass

    level1.aoa_mask = aoa_mask
    level1.vz_mask = vz_mask

    return aoa_mask, vz_mask, aoa_store


def _dy_dx(y_arr, x_arr):
    out = np.zeros(np.shape(y_arr))
    for x1, x2, y1, y2, i in \
            zip(x_arr[0:-1], x_arr[1:], y_arr[0:-1], y_arr[1:], [i+1 for i in (range(y_arr.shape[0]-1))]):
        out[i] = (y2 - y1) / (x2 - x1)
    return out


def _mag(vector):
    return np.sqrt(np.sum(vector**2))


def detect_cloud_limits(level1, profile, prof_num=1, prom=10, diff_lim=0.00002, detect_type="cursor"):

    if "CYISUAData" in str(type(level1)):
        conc = level1.mass_concentration1
    else:
        conc = level1.mass_concentration
    if profile == "up":
        mask = level1.up_profile_mask
    elif profile == "down":
        mask = level1.down_profile_mask
    else:
        raise ValueError("ERROR: Invalid profile type, up or down are the options")
    t = level1.time
    conc = conc[np.where(mask[:, prof_num-1] == 1)]

    conc_diff = _dy_dx(conc, t)
    norm_conc_diff = conc_diff.astype(float) - float(conc_diff[0])
    limits = []

    if detect_type == "peaks":
        norm_conc_diff = norm_conc_diff[norm_conc_diff > -800]
        p_peaks, _ = find_peaks(np.squeeze(norm_conc_diff), prominence=prom, distance=10)
        n_peaks, _ = find_peaks(np.squeeze(norm_conc_diff) * -1, prominence=prom, distance=10)

    elif detect_type == "cursor":
        head_cursor = 0
        for conc_diff_val in norm_conc_diff:
            if abs(conc_diff_val) >= abs(diff_lim):
                limits.append(head_cursor)
                break
            head_cursor += 1
        tail_cursor = int(np.shape(norm_conc_diff)[0])
        for conc_diff_val in np.flip(np.squeeze(norm_conc_diff)):
            if abs(conc_diff_val) >= abs(diff_lim):
                limits.append(tail_cursor)
                break
            tail_cursor -= 1

    else:
        raise ValueError("ERROR: Detect type is invalid")

    if len(limits) != 2:
        raise ValueError("ERROR: no limits detected, diff_lim is too high")

    return limits


def adjust_airspeed_mtof(level1, profile, window=15, prof_num=1, prom=10, diff_lim=0.00002, detect_type="cursor"):
    try:
        _ = level1.adjusted_airspeed
    except(NameError, AttributeError):
        raise AttributeError("ERROR: No adjusted_airspeed, data may not be FW")

    if "CYISUAData" in str(type(level1)):
        mtof = np.true_divide(80, np.true_divide(level1.m_tof1[:, 0], 3))
    else:
        mtof = np.true_divide(80, np.true_divide(level1.m_tof[:, 0], 3))
    if profile == "up":
        mask = level1.up_profile_mask
    elif profile == "down":
        mask = level1.down_profile_mask
    else:
        raise ValueError("ERROR: Invalid profile type, up or down are the options")

    x_mask = np.array(range(mtof.shape[0]))[np.isfinite(mtof)]
    y_mask = mtof[np.isfinite(mtof)]
    mtof_interp = interp1d(x_mask, y_mask, fill_value="extrapolate")
    mtof = mtof_interp(range(mtof.shape[0]))
    mtof = common.auto_conv_filter(window, mtof)

    airspeed = np.true_divide(level1.vz_cms, 100.0)
    mtof = mtof[np.where(mask[:, prof_num-1] == 1)]
    airspeed = airspeed[np.where(mask[:, prof_num-1] == 1)]
    cloud_limits = detect_cloud_limits(level1, profile, prom=prom, diff_lim=diff_lim, detect_type=detect_type)

    cal1 = airspeed[cloud_limits[0]]
    cal2 = airspeed[cloud_limits[1]]

    off1 = cal1 - mtof[cloud_limits[0]]
    off2 = cal2 - mtof[cloud_limits[1]]

    m = np.true_divide((off2 - off1), (cloud_limits[1] - cloud_limits[0]))
    c = off1

    adj_speed = np.zeros(np.shape(airspeed))
    index = 0
    for s, i in zip(airspeed, range(airspeed.shape[0])):
        if cloud_limits[0] <= i <= cloud_limits[1]:
            adj_speed[i] = mtof[i] + (m*index + c)
            index = index + 1
        else:
            adj_speed[i] = s

    level1.adjusted_airspeed = adj_speed

    return adj_speed
