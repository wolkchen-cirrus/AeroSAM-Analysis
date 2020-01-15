import numpy as np
import cPickle as Pickle


def import_level1(level1_path):
    with open(level1_path) as level1_file:
        level1_data = Pickle.load(level1_file)
    return level1_data


def fetch_row(altitude=None, time=None, level1_data=None):

    a = str(type(level1_data))

    if "SUAData" in str(type(level1_data)):
        try:
            key_col = level1_data.alt
            row_value = altitude
        except AttributeError:
            raise AttributeError("ERROR: level1_data object problem")
    elif ("StaticCASData" in str(type(level1_data))) or ("StaticFSSPData" in str(type(level1_data))):
        try:
            key_col = level1_data.time
            row_value = time
        except AttributeError:
            raise AttributeError("ERROR: level1_data object problem")
    else:
        raise ValueError("ERROR: Unrecognised data object")
    if not isinstance(key_col, np.ndarray):
        try:
            key_col = np.asarray(key_col)
        except(ValueError, TypeError):
            raise TypeError("ERROR: Incompatible data type")
    diff_col = abs(key_col - row_value)
    min_diff = np.amin(diff_col)
    min_diff_index = np.where(diff_col == min_diff)

    return min_diff_index
