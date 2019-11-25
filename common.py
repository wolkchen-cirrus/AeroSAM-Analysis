import os
import numpy as np
import calendar


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
