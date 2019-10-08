"""
This script contains all the functions for converting level 0 SUA data into level 1 SUA data. Each function
requires a "SUAData" object from the "importer.py script
"""

import numpy as np
from scipy.signal import find_peaks
import common
from os import listdir
from os import name as osname


def split_by_pressure(sua_data):
    """
    This function will create two "profile masks" for upwards and downwards data. This is an ndarray of 1s and 0s, the
    columns of which are the same size as the columns of data imported from the csv files. Multiply a column from the
    mask by a column of data to obtain a profile. Pressure is used as the variable to determine when the profiles are.
    Different columns of the mask mean different profiles in one csv file.
    :param sua_data: The SUA data object from "importer.py".
    :return: Nothing, all assignments effect input object.
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


def assign_ucass_lut(sua_data, material="Water", path=None):
    """
    This function chooses the correct look up table (LUT) for the UCASS used in the SUA data. The LUT is chosen firstly
    based on tags similarity, and second based on date (closest date to that of the SUA data object). The LUT is stored
    as a dict object. The matching tags for the LUT are specified in the LUT's filename.
    :param sua_data: The SUA data object from "importer.py".
    :param material: The material used in the LUT generation (Water by default).
    :param path: If specified, the function will assign a user defined path rather than inferring from tags.
    :return: Nothing, all assignments effect input object.
    """

    # Ensuring there are no problems with the SUA data class and importing the tags and date.
    try:
        tags = sua_data.tags
        data_date = sua_data.datetime
    except NameError:
        raise NameError("ERROR: Problem with SUA data object")
    if sua_data.ucass_lut_aerosol is not None:
        print "WARNING: Overwriting existing UCASS LUT"
    elif sua_data.ucass_lut_droplet is not None:
        print "WARNING: Overwriting existing UCASS LUT"

    date = int(data_date.split()[0].replace("-", ""))   # Format date into computer readable format
    tags.append(material)                               # Add the material as a tag

    # If the Uncalibrated tag exists, there will be no relevant LUT.
    if "Uncalibrated" in tags:
        print "WARNING: UCASS is uncalibrated, no LUT could be assigned"
        return

    # If the path is specified, assign this directly and do not search tags.
    elif path:
        if "Aerosol" in tags:
            print "INFO: Aerosol sonde detected, LUT at user specified path"
            sua_data.ucass_lut_aerosol = common.file_to_dict(path)
        if "Droplet" in tags:
            print "INFO: Droplet sonde detected, LUT at user specified path"
            sua_data.ucass_lut_droplet = common.file_to_dict(path)
        return

    # The code for tag searching and similarity computation.
    else:
        # First, compute how many tags each LUT file shares in common with the SUA data object.
        lut_dir_path = common.read_setting("lut_dir_path")  # LUT directory defined in settings
        lut_files = listdir(lut_dir_path)                   # Get files in directory
        similarity = []                                     # Pre-assign similarity list
        index = 0                                           # Index of file, order is constant
        for lut in lut_files:                               # Cycle through LUT files
            similarity.append(0)                            # Start with no similar tags, then add one if detected
            lut_tags = lut.split('_')                       # Tags in file name delimited with _
            for lut_tag in lut_tags:                        # Cycle through tags in LUT file name
                if lut_tag in tags:                         # Check if each tag is in the SUA data tags
                    similarity[index] += 1                  # If yes, increase similarity by 1
            index += 1
        lut_index = max(similarity)                         # Find max similarity

        # In the case that equal max similarities are detected, use date as the deciding factor. The date on the LUT
        # should be chosen to be the most similar to the date on the SUA data object.
        if similarity.count(lut_index) > 1:     # Check if there is multiple max similarities
            date_list = []                      # Pre assign list of dates

            # Fill up the list of dates using the last tag in the LUT filename (minus extension).
            for lut in lut_files:
                lut_date = int(lut.split('_')[-1].replace(".LUT", ""))
                date_list.append(lut_date)

            # Compute which date is the most similar to the SUA data date.
            date_list = [abs(x - date) for x in date_list]  # Subtract SUA data date from the LUT date
            mask = []                                       # Mask out all the elements with low similarity
            for i in similarity:
                if i == lut_index:
                    mask.append(1)
                else:
                    mask.append(0)
            date_list = [mask[i] * date_list[i] for i in range(len(mask))]  # Apply the mask to the list of dates
            chosen_lut_date = min(date_list)                                # Find minimum
            lut = date_list.index(chosen_lut_date)                          # Find location of minimum
        else:
            lut = similarity.index(lut_index)                               # LUT index if there is only one LUT

        # Assign the correct LUT file to the SUA data object
        lut_file = lut_files[lut]                       # Get file name
        lut_path = ""
        if osname == 'nt':                                # If windows
            lut_path = lut_dir_path + "\\" + lut_file
        elif osname == 'posix':                           # If Linux
            lut_path = lut_dir_path + "/" + lut_file

        if "Aerosol" in tags:
            print "INFO: Aerosol sonde detected, LUT chosen was: %s" % lut_file
            sua_data.ucass_lut_aerosol = common.file_to_dict(lut_path)
        if "Droplet" in tags:
            print "INFO: Droplet sonde detected, LUT chosen was: %s" % lut_file
            sua_data.ucass_lut_droplet = common.file_to_dict(lut_path)

    return
