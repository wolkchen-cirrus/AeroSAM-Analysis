"""
level0to1: This script contains all the functions for converting level 0 SUA data into level 1 SUA data. Each function
requires a "SUAData" object from the "importer.py script
"""

import numpy as np
from scipy.signal import find_peaks
import common


def split_by_vz(sua_data):
    vz_lim = float(common.read_setting("vz_lim"))
    try:
        vz_cms = sua_data.vz_cms
    except NameError:
        raise NameError("ERROR: Problem with SUA data object")

    return
