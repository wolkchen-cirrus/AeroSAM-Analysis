import os
import numpy as np


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


class ColumnProperty(object):
    def __init__(self, name):
        self.name = "_" + name

    def __get__(self, obj, cls=None):
        return getattr(obj, self.name)

    def __set__(self, obj, value):
        if isinstance(value, str):
            try:
                float(value)
            except TypeError:
                raise TypeError("ERROR: Invalid Type")
        if isinstance(value, list):
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


class SplitCounter(object):
    def __init__(self, name):
        self.name = "_" + name

    def __get__(self, obj, cls=None):
        return getattr(obj, self.name)

    def __set__(self, obj, value):

        try:
            value = int(value)
        except TypeError:
            raise TypeError("ERROR: Invalid Type")

        if getattr(obj, self.name) is None:
            print "INFO: Setting fist value"
            setattr(obj, self.name, np.array([value]))
        else:

            if value:
                current = getattr(obj, self.name)
                flat_current = current*0
                i = 0
                for row in current:
                    try:
                        flat_current[i] = row/row
                    except (ZeroDivisionError, RuntimeWarning):
                        continue
                    i += 1

                flat_current[-1] = value

                diff = current*0
                i = 0
                for row in flat_current:
                    try:
                        diff[i] = flat_current[i+1] - row
                        if diff[i] < 0:
                            diff[i] = 0
                        i += 1
                    except IndexError:
                        pass

                placeholder = []
                for row in diff:
                    if row:
                        placeholder.append(row)

                corrected_value = int(len(placeholder))

            else:
                corrected_value = value

            corrected_value = np.array([corrected_value])
            new_value = np.vstack([getattr(obj, self.name), corrected_value])
            filtered_value = self._zero_filter(new_value, 3, 7)
            setattr(obj, self.name, filtered_value)

    @classmethod
    def _zero_filter(cls, data, window_size, threshold):
        new_data = data
        i = 0
        buf = []
        for row in data:
            try:
                if (row == data[i+1]) and (row is not 0):
                    buf.append(i)
                i += 1
            except IndexError:
                pass
        for i in buf:
            window = buf[i:(i+window_size)]
            normal_window = [x - window[0] for x in window]
            if sum(normal_window) > threshold:
                for index in window:
                    new_data[index] = 0
        return new_data


def line_nums(file_name):
    with open(file_name) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
