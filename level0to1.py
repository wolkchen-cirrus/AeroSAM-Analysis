"""
This script contains all the functions for converting level 0 SUA data into level 1 SUA data. Each function
requires a "SUAData" object from the "importer.py script
"""

import numpy as np
from scipy.signal import find_peaks
import common


def split_by_pressure(sua_data):
    """
    This function will create two "profile masks" for upwards and downwards data. This is an ndarray of 1s and 0s, the
    columns of which are the same size as the columns of data imported from the csv files. Multiply a column from the
    mask by a column of data to obtain a profile. Pressure is used as the variable to determine when the profiles are.
    Different columns of the mask mean different profiles in one csv file.
    :param sua_data: The SUA data object from "importer.py"
    :return: Nothing
    """

    # Ensuring there are no problems with the SUA data class and importing the pressure.
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

    # Using SciPy to find the peaks in the pressure, used as an indicator of how many profiles are in the data set.
    norm_press_hpa = press_hpa.astype(float) - float(press_hpa[0])          # Normalizing pressure
    p_peaks, _ = find_peaks(np.squeeze(norm_press_hpa), prominence=1)       # Positive peaks in data
    n_peaks, _ = find_peaks(np.squeeze(norm_press_hpa) * -1, prominence=1)  # Negative peaks in data
    num_peaks = int(np.shape(n_peaks)[0])                                   # Number of profiles

    # Removing the waiting time, which manifests as a head and tail in the data.
    press_lim = float(common.read_setting("press_lim"))                 # Pressure limit for motion
    head_cursor = 0                                                     # Cursor for data head
    for pressure in norm_press_hpa:                                     # Cycling through till the limit is reached
        if abs(pressure) >= abs(press_lim):
            p_peaks = np.hstack((np.array([head_cursor]), p_peaks))     # Add cursor index to p_peaks
            break
        head_cursor += 1
    tail_cursor = int(np.shape(norm_press_hpa)[0])                      # Tail cursor starts at the end
    for pressure in np.flip(np.squeeze(norm_press_hpa)):                # Count from the back of array till limit
        if abs(pressure) >= abs(press_lim):
            p_peaks = np.hstack((p_peaks, np.array([tail_cursor])))     # append index to end of p_peaks
            break
        tail_cursor -= 1                                                # subtract from tail cursor each loop

    # Check here to ensure the number of positive peaks is one more than the number of negative peaks (what goes up,
    # must come down!)
    if not (int(np.shape(n_peaks)[0]) + 1) == int(np.shape(p_peaks)[0]):
        raise ValueError("ERROR: Problem detecting peaks")

    # Creating the profile masks from the detected peak indices.
    up_profile_store = np.zeros((int(np.shape(norm_press_hpa)[0]), num_peaks))
    down_profile_store = np.zeros((int(np.shape(norm_press_hpa)[0]), num_peaks))
    for i in range(num_peaks):
        up_profile_store[p_peaks[i]:n_peaks[i], i] = 1
        down_profile_store[n_peaks[i]:p_peaks[i+1], i] = 1

    # Assigning the profile masks to the SUA data properties.
    try:
        sua_data.up_profile_mask = up_profile_store
        sua_data.down_profile_mask = down_profile_store
    except NameError:
        raise NameError("ERROR: Problem with SUA data object")

    return
