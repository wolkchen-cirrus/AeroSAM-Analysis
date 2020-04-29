import numpy as np
import cPickle as Pickle
from AirborneParticleAnalysis import common
from scipy import stats


def import_level1(level1_path):
    with open(level1_path) as level1_file:
        level1_data = Pickle.load(level1_file)
    return level1_data


def fetch_row(altitude=None, time=None, level1_data=None):

    if "SUAData" in str(type(level1_data)):
        try:
            key_col = level1_data.alt
            if not isinstance(key_col, np.ndarray):
                try:
                    key_col = np.asarray(key_col)
                except(ValueError, TypeError):
                    raise TypeError("ERROR: Incompatible data type")
            row_value = altitude
            prof_mask = level1_data.up_profile_mask
            (r, cols) = prof_mask.shape
            rows = []
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


def fetch_row_tolerance(altitude=None, time=None, level1_data=None):

    if ("SUAData" in str(type(level1_data))) or ("CYISUAData" in str(type(level1_data))):
        try:
            key_col = level1_data.alt
            tol = float(common.read_setting("height_mean_tolerance_metres"))*1000.0
            if not isinstance(key_col, np.ndarray):
                try:
                    key_col = np.asarray(key_col)
                except(ValueError, TypeError):
                    raise TypeError("ERROR: Incompatible data type")
            row_value = altitude
            prof_mask = level1_data.up_profile_mask
            (r, cols) = prof_mask.shape
            rows = []
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

                rows.append(buf)
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

            rows.append(buf)
        except AttributeError:
            raise AttributeError("ERROR: level1_data object problem")

    else:
        raise ValueError("ERROR: Unrecognised data object")

    return rows


def mean_dn_dlogdp(level1_data, rows, ucass_number=1):

    if isinstance(rows[0], list):
        raise ValueError("ERROR: Pass only one profile into function")

    if "CYISUAData" in str(type(level1_data)):

        if ucass_number == 1:
            dn_dlogdp = level1_data.dn_dlogdp1
        elif ucass_number == 2:
            dn_dlogdp = level1_data.dn_dlogdp2
        else:
            raise ValueError("ERROR: For CYI, only 2 UCASS' are configured")

        data_arr = np.zeros([len(rows), len(dn_dlogdp[rows[0]])])
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
