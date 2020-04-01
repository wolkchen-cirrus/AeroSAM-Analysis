import os
import numpy as np
import calendar
import math


def read_setting(setting):

    module_path = os.path.dirname(os.path.realpath(__file__))
    settings_path = module_path + "/settings.txt"
    if os.path.exists(settings_path):
        pass
    else:
        raise ValueError("ERROR: Settings path does not exist")
    try:
        settings_file = open(settings_path, 'r')
    except ValueError:
        raise ValueError("ERROR: Settings file could not be opened")
    set_line = []

    if settings_file.mode == 'r':
        settings_string = settings_file.read()
        settings_list = settings_string.split('\n')
        for i in settings_list:
            if setting in i:
                set_line = i
    else:
        settings_file.close()
        raise RuntimeError("ERROR: File not opened")

    if not set_line:
        settings_file.close()
        raise RuntimeError("ERROR: Setting not found")
    else:
        set_line_list = set_line.split(' = ')
        settings_file.close()
        return set_line_list[1]


class AddedColumn(object):
    def __init__(self, name):
        self.name = "_" + name

    def __get__(self, obj, cls=None):
        return getattr(obj, self.name)

    def __set__(self, obj, value):
        try:
            value = value.astype(float)
        except TypeError:
            raise TypeError("ERROR: Invalid Type")

        line_num = getattr(obj, "num_lines")
        size = value.shape

        if size[0] != line_num:
            raise ValueError("ERROR: Array is not the same column length")

        setattr(obj, self.name, value)


class ColumnProperty(object):
    def __init__(self, name):
        self.name = "_" + name

    def __get__(self, obj, cls=None):
        return getattr(obj, self.name)

    def __set__(self, obj, value):
        if isinstance(value, str):
            try:
                value = float(value)
            except TypeError:
                raise TypeError("ERROR: Invalid Type")
        if isinstance(value, list):
            value = map(int, value)
            if getattr(obj, self.name) is None:
                print "INFO: Setting fist value"
                setattr(obj, self.name, np.array(value))
            else:
                value = np.array(value)
                new_value = np.vstack([getattr(obj, self.name), value])
                setattr(obj, self.name, new_value)
        else:
            if getattr(obj, self.name) is None:
                print "INFO: Setting first value"
                setattr(obj, self.name, np.array([value]))
            else:
                value = np.array([value])
                new_value = np.vstack([getattr(obj, self.name), value])
                setattr(obj, self.name, new_value)


def line_nums(file_name):
    with open(file_name) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def file_to_dict(path):
    d = {}
    with open(path) as f:
        for line in f:
            (key, val) = line.split(',')
            try:
                d[int(key)] = float(val.replace('\n', ''))
            except ValueError:
                try:
                    d[key] = float(val.replace('\n', ''))
                except ValueError:
                    pass
    return d


def utc_to_epoch(timestamp):
    epoch = calendar.timegm(timestamp.utctimetuple())
    return epoch


def make_file(file_path, extension, base_name=None):
    if base_name:
        name = base_name
        path = file_path
    else:
        name = file_path.split("\\")[-1].split(".")[0].split("_")[0:-1]
        path_l = file_path.split("\\")
        del path_l[-1]
        path = "\\".join(path_l)
    name = "_".join(name)
    name += "_00"
    for i in range(100):
        path_name = path
        name_l = list(name)
        name_l[-1 - 1] = str(int(i / 10))
        name_l[-1] = str(int(i % 10))
        name = "".join(name_l)
        path_name += '\\'
        path_name += name
        path_name += extension
        if os.path.exists(path_name) is False:
            return path_name


def seconds_to_timestamp(seconds):
    hours = int(math.floor(seconds / 3600.0))
    mins = int(math.floor((seconds % 3600.0) / 3600.0 * 60))
    secs = int((seconds % 3600.0) % 60)
    hhmmss = str(hours) + ":" + str(mins) + ":" + str(secs)
    return hhmmss


def hhmmss_to_sec(hhmmss):
    return int(hhmmss[0:2])*3600+int(hhmmss[2:4])*60+int(hhmmss[4:6])


def cm_to_inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(k/inch for k in tupl[0])
    else:
        return tuple(k/inch for k in tupl)


def sync_data(t1, t2, t_list_hr, d_list_hr, sync_type="mean"):

    if not isinstance(t_list_hr, float):
        try:
            t_list_hr = float(t_list_hr)
        except (ValueError, TypeError):
            raise TypeError("ERROR: t_list_hr must be floats")
    elif not isinstance(d_list_hr, float):
        try:
            d_list_hr = float(d_list_hr)
        except (ValueError, TypeError):
            raise TypeError("ERROR: d_list_hr must be floats")
    elif not isinstance(t1, float):
        try:
            t1 = float(t1)
        except (ValueError, TypeError):
            raise TypeError("ERROR: t1 must be float")
    elif not isinstance(t2, float):
        try:
            t2 = float(t2)
        except (ValueError, TypeError):
            raise TypeError("ERROR: t2 must be float")

    index = 0
    t1_index = None
    t2_index = None
    for t1_hr, t2_hr in zip(t_list_hr[0:-1], t_list_hr[1:]):

        if t1_hr >= t1 > t2_hr:
            t1_index = index
        elif t1_hr >= t2 > t2_hr:
            t2_index = index

        if (t1_index is not None) and (t2_index is not None):
            break

        index += 1

    if (t1_index is None) or (t2_index is None):
        raise ValueError("ERROR: Index not assigned, check function inputs")

    d_range = d_list_hr[t1_index:t2_index]

    if sync_type == "mean":
        return sum(d_range)/float(len(d_range))
    elif sync_type == "sum":
        return sum(d_range)
    else:
        raise ValueError("ERROR: Valid sync types are \"mean\" and \"sum\"")
