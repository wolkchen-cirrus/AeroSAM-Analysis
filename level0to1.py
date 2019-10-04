"""
level0to1: This script contains all the functions for converting level 0 SUA data into level 1 SUA data. Each function
requires a "SUAData" object from the "importer.py script
"""

import numpy as np
from scipy.signal import find_peaks
import common


def split_by_pressure(sua_data):
    try:
        press_hpa = sua_data.press_hpa
    except NameError:
        raise NameError("ERROR: Problem with SUA data object")
    if press_hpa is None:
        raise ValueError("ERROR: SUA data is not level 0")
    if sua_data.up_profile_mask is not None:
        print "WARNING: Overwriting existing profile analysis"
    elif sua_data.down_profile_mask is not None:
        print "WARNING: Overwriting existing profile analysis"

    norm_press_hpa = press_hpa.astype(float) - float(press_hpa[0])
    p_peaks, _ = find_peaks(np.squeeze(norm_press_hpa), prominence=1)
    n_peaks, _ = find_peaks(np.squeeze(norm_press_hpa) * -1, prominence=1)
    num_peaks = int(np.shape(n_peaks)[0])

    press_lim = float(common.read_setting("press_lim"))
    head_cursor = 0
    for pressure in norm_press_hpa:
        if abs(pressure) >= abs(press_lim):
            p_peaks = np.hstack((np.array([head_cursor]), p_peaks))
            break
        head_cursor += 1
    tail_cursor = int(np.shape(norm_press_hpa)[0])
    for pressure in np.flip(np.squeeze(norm_press_hpa)):
        if abs(pressure) >= abs(press_lim):
            p_peaks = np.hstack((p_peaks, np.array([tail_cursor])))
            break
        tail_cursor -= 1

    if not (int(np.shape(n_peaks)[0]) + 1) == int(np.shape(p_peaks)[0]):
        raise ValueError("ERROR: Problem detecting peaks")

    up_profile_store = np.zeros((int(np.shape(norm_press_hpa)[0]), num_peaks))
    down_profile_store = np.zeros((int(np.shape(norm_press_hpa)[0]), num_peaks))
    for i in range(num_peaks):
        up_profile_store[p_peaks[i]:n_peaks[i], i] = 1
        down_profile_store[n_peaks[i]:p_peaks[i+1], i] = 1

    try:
        sua_data.up_profile_mask = up_profile_store
        sua_data.down_profile_mask = down_profile_store
    except NameError:
        raise NameError("ERROR: Problem with SUA data object")

    return
