"""
This script contains all the functions for converting level 0 SUA data into level 1 SUA data. Each function
requires a "SUAData" object from the "importer.py script.
"""

import numpy as np
from scipy.signal import find_peaks
from AirborneParticleAnalysis import common
from os import listdir
from os import name as osname
from os import path as ospath
import cPickle as Pickle
import warnings


def export_level1(level1_data):
    """
    This function pickles the level 1 data object and saves it to a file for easier imports in the future. It also
    checks to ensure data is at level 1 before exporting.
    :param level1_data: The level 1 data object
    :return: Nothing, all assignments effect input object.
    """

    # Ensuring there are no problems with the SUA data class.
    try:
        current_path = level1_data.path
        level1_data.check_level()                   # Run the level check
        data_level = level1_data.level_indicator    # Get data level
    except NameError:
        raise NameError("ERROR: Problem with SUA data object")

    # Makes sure the data is level 1 before continuing the export
    if data_level != 1:
        raise ValueError("ERROR: Specified data must be level 1 or above")

    # Getting the level 1 directory path
    path_list = current_path.split("\\")
    path_list[-1-1] = path_list[-1-1].replace("0", "1")
    level1_data.path = "\\".join(path_list)
    path = level1_data.path
    file_name = common.make_file(path, ".pdat")         # Make the level 1 file name
    with open(file_name, "w+") as out_file:
        Pickle.dump(level1_data, out_file)              # Saving the pickled object to the file
    return


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
        if "CYISUAData" in str(type(sua_data)):
            press_hpa = np.asarray(press_hpa)
            press_hpa = press_hpa[press_hpa != 0]
    except NameError:
        raise NameError("ERROR: Problem with SUA data object")
    if press_hpa is None:
        raise ValueError("ERROR: SUA data is not level 0")
    if sua_data.up_profile_mask is not None:
        warnings.warn("WARNING: Overwriting existing profile analysis")
    elif sua_data.down_profile_mask is not None:
        warnings.warn("WARNING: Overwriting existing profile analysis")

    if "FMISUAData" in str(type(sua_data)):
        prom = 5
    else:
        prom = 1

    # Using SciPy to find the peaks in the pressure, used as an indicator of how many profiles are in the data set.
    norm_press_hpa = press_hpa.astype(float) - float(press_hpa[0])                              # Normalizing pressure
    norm_press_hpa = norm_press_hpa[norm_press_hpa > -800]
    p_peaks, _ = find_peaks(np.squeeze(norm_press_hpa), prominence=prom, distance=10)           # Positive peaks in data
    n_peaks, _ = find_peaks(np.squeeze(norm_press_hpa) * -1, prominence=prom, distance=10)      # Negative peaks in data
    num_peaks = int(np.shape(n_peaks)[0])                                                       # Number of profiles

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
        name_arr = None
        if "CYISUAData" in str(type(sua_data)):
            ucass_name1 = sua_data.ucass_name1
            ucass_name2 = sua_data.ucass_name2
            name_arr = [ucass_name1, ucass_name2]
    except NameError:
        raise NameError("ERROR: Problem with SUA data object")

    date = int(data_date.split()[0].replace("-", ""))   # Format date into computer readable format
    tags.append(material)                               # Add the material as a tag
    if "CYISUAData" in str(type(sua_data)):
        tags.append("xx")

    # If the Uncalibrated tag exists, there will be no relevant LUT.
    if "Uncalibrated" in tags:
        warnings.warn("WARNING: UCASS is uncalibrated, no LUT could be assigned")
        return

    # If the path is specified, assign this directly and do not search tags.
    elif path:
        if ("SUAData" in str(type(sua_data))) and ("CYI" not in str(type(sua_data))):
            if "Aerosol" in tags:
                print "INFO: Aerosol sonde detected, LUT at user specified path"
                sua_data.ucass_lut_aerosol = common.file_to_dict(path)
            if "Droplet" in tags:
                print "INFO: Droplet sonde detected, LUT at user specified path"
                sua_data.ucass_lut_droplet = common.file_to_dict(path)
        elif "CYISUAData" in str(type(sua_data)):
            raise NotImplementedError("ERROR: Path specification not supported for CYISUAData")
        return

    # The code for tag searching and similarity computation.
    else:
        if "CYISUAData" in str(type(sua_data)):
            ucass_amount = 2
        else:
            ucass_amount = 1

        # First, compute how many tags each LUT file shares in common with the SUA data object.
        similarity_arr = []
        lut_dir_path = common.read_setting("lut_dir_path")      # LUT directory defined in settings
        lut_files = listdir(lut_dir_path)                       # Get files in directory
        lut_index = []
        for i in range(ucass_amount):
            similarity = []                                     # Pre-assign similarity list
            index = 0                                           # Index of file, order is constant
            if "CYISUAData" in str(type(sua_data)):
                tags[-1] = name_arr[i]
            for lut in lut_files:                               # Cycle through LUT files
                similarity.append(0)                            # Start with no similar tags, then add one if detected
                lut_tags = lut.split('_')                       # Tags in file name delimited with _
                for lut_tag in lut_tags:                        # Cycle through tags in LUT file name
                    if lut_tag in tags:                         # Check if each tag is in the SUA data tags
                        similarity[index] += 1                  # If yes, increase similarity by 1
                index += 1
            similarity_arr.append(similarity)
            lut_index.append(max(similarity))                   # Find max similarity

        # In the case that equal max similarities are detected, use date as the deciding factor. The date on the LUT
        # should be chosen to be the most similar to the date on the SUA data object.
        do_date_loop = None
        for i, n in zip(lut_index, range(len(lut_index))):
            if similarity_arr[n].count(i) > 1:                     # Check if there is multiple max similarities
                do_date_loop = 1

        if do_date_loop:
            lut = []
            for j, n in zip(lut_index, range(len(lut_index))):
                index = 0
                lut_candidate_index = []
                for i in similarity_arr[n]:
                    if i == j:
                        lut_candidate_index.append(index)
                    index += 1
                index = 0
                new_luts = []
                for i in lut_files:
                    if index in lut_candidate_index:
                        new_luts.append(i)
                    index += 1

                date_list = []                      # Pre assign list of dates

                # Fill up the list of dates using the last tag in the LUT filename (minus extension).
                for i in new_luts:
                    lut_date = int(i.split('_')[-1].replace(".LUT", ""))
                    date_list.append(lut_date)

                # Compute which date is the most similar to the SUA data date.
                date_list = [abs(x - date) for x in date_list]              # Subtract SUA data date from the LUT date
                chosen_lut_date = min(date_list)                            # Find minimum
                new_lut_index = date_list.index(chosen_lut_date)            # Find location of minimum
                chosen_lut = new_luts[new_lut_index]
                lut.append(lut_files.index(chosen_lut))
        else:
            lut = []
            for i, n in zip(lut_index, range(len(lut_index))):
                lut.append(similarity_arr[n].index(i))                             # LUT index if there is only one LUT

        # Assign the correct LUT file to the SUA data object
        lut_paths = []
        for i in lut:
            lut_file = lut_files[i]                             # Get file name
            lut_path = ""
            if osname == 'nt':                                  # If windows
                lut_path = lut_dir_path + "\\" + lut_file
            elif osname == 'posix':                             # If Linux
                lut_path = lut_dir_path + "/" + lut_file
            lut_paths.append(lut_path)

        if "CYISUAData" in str(type(sua_data)):
            if str(lut_paths[0].split("\\")[-1].split("_")[0]) in sua_data.ucass_name1:
                sua_data.ucass_lut1 = common.file_to_dict(lut_paths[0])
                if "Droplet" in str(lut_paths[0].split("\\")[-1]):
                    sua_data.ucass_gain1 = "Droplet"
                elif "Aerosol" in str(lut_paths[0].split("\\")[-1]):
                    sua_data.ucass_gain1 = "Aerosol"
                else:
                    raise ValueError("ERROR: No gain in LUT")
            else:
                raise ValueError("ERROR: Invalid UCASS LUT Assigned")
            if str(lut_paths[1].split("\\")[-1].split("_")[0]) in sua_data.ucass_name2:
                sua_data.ucass_lut2 = common.file_to_dict(lut_paths[1])
                if "Droplet" in str(lut_paths[0].split("\\")[-1]):
                    sua_data.ucass_gain2 = "Droplet"
                elif "Aerosol" in str(lut_paths[0].split("\\")[-1]):
                    sua_data.ucass_gain2 = "Aerosol"
                else:
                    raise ValueError("ERROR: No gain in LUT")
            else:
                raise ValueError("ERROR: Invalid UCASS LUT Assigned")
        else:
            for i in lut_paths:
                if "Aerosol" in tags:
                    print "INFO: Aerosol sonde detected, LUT chosen was: %s" % i
                    sua_data.ucass_lut_aerosol = common.file_to_dict(i)
                elif "Droplet" in tags:
                    print "INFO: Droplet sonde detected, LUT chosen was: %s" % i
                    sua_data.ucass_lut_droplet = common.file_to_dict(i)

    return


def bin_centre_dp_um(sua_data, ignore_b1=False, centre_type="Geometric"):
    """
    This function will compute the bin centres according to some user-specified options. This must be done before all
    particle counter data analysis.
    :param sua_data: The data object, note the code will only accept certain types
    :param ignore_b1: False by default, will effectively ignore the first bin (common in some analyses)
    :param centre_type: The averaging type used to compute centring, will accept Geometric or Arithmetic
    :return: Nothing, all assignments effect input object.
    """

    # Ensuring there are no problems with the SUA data class and importing.
    tags = None
    ubs = None
    ubs1 = None
    ubs2 = None
    try:
        if "CYISUAData" in str(type(sua_data)):
            ubs1 = sua_data.bins1
            ubs2 = sua_data.bins2
        else:
            ubs = sua_data.bins
            tags = sua_data.tags
    except NameError:
        raise NameError("ERROR: Problem with SUA data object")

    # Checking if the input object type is SUAData, analysis will depend on object type. The first stage on computation
    # is converting the bin boundaries from 12 bit ADC to a diameter in um. Note the CAS data does not need this since
    # it is already listed at level 0.
    if ("SUAData" in str(type(sua_data))) and ("CYI" not in str(type(sua_data))):
        ubs_um_list = []

        # If no UCASS LUT is assigned, this function must be run first
        if (sua_data.ucass_lut_aerosol is None) and (sua_data.ucass_lut_droplet is None):
            print("INFO: Assigning LUT to UCASS")
            assign_ucass_lut(sua_data)

        # Check if the UCASS is in droplet or aerosol mode and get the lookup table. The 'ucass_type' variable is used
        # because multiple properties depend on the UCASS gain mode.
        if "Droplet" in tags:
            lut = sua_data.ucass_lut_droplet
            ucass_type = 0
        elif "Aerosol" in tags:
            lut = sua_data.ucass_lut_aerosol
            ucass_type = 1
        else:
            raise ValueError("ERROR: Gain not specified in tags, cannot compute bin centres")

        # Convert instrument response (12 bit ADC) into a size using the lookup table, this is a dict object where the
        # key is ADC and the choresponding value is diameter in um.
        ubs_um = []
        for i in ubs:
            ubs_um.append(lut[int(i)])

        # Assign the first bin lower boundary if the ignore_b1 flag is False.
        if ignore_b1 is False:
            if ucass_type == 0:
                b1_lb = 1.0             # For low gain (Droplet)
            elif ucass_type == 1:
                b1_lb = 0.3             # For high gain (Aerosol)
            else:
                raise ValueError("ERROR: Invalid UCASS gain mode")
            ubs_um.insert(0, b1_lb)     # Insert value at start of list
        elif ignore_b1 is True:
            pass
        else:
            raise ValueError("ERROR: \'ignore_b1\' is not boolean")

        ubs_um = [float(i) for i in ubs_um]     # Convert to floats for analysis
        sua_data.bin_bounds_dp_um = ubs_um      # Assign bin bounds property (legacy)
        ubs_um_list.append(ubs_um)

    # Check if the input object is Static CAS data, in which case the above need not be performed.
    elif ("StaticCASData" in str(type(sua_data))) or ("StaticFSSPData" in str(type(sua_data))):
        ubs_um_list = []
        if ignore_b1:
            del ubs[0]                          # Delete first bin if ignore_b1 flag, this may cause index errors
        ubs_um = ubs
        ubs_um = [float(i) for i in ubs_um]     # Convert to float
        ubs_um_list.append(ubs_um)

    elif "CYISUAData" in str(type(sua_data)):
        # If no UCASS LUT is assigned, this function must be run first
        if (sua_data.ucass_lut1 is None) or (sua_data.ucass_lut2 is None):
            print("INFO: Assigning LUT to UCASS")
            assign_ucass_lut(sua_data)

        lut_arr = [sua_data.ucass_lut1, sua_data.ucass_lut2]
        gain_arr = [sua_data.ucass_gain1, sua_data.ucass_gain2]
        ubs_arr = [ubs1, ubs2]

        # Convert instrument response (12 bit ADC) into a size using the lookup table, this is a dict object where the
        # key is ADC and the choresponding value is diameter in um.
        ubs_um_list = []
        for lut, gain, ubs in zip(lut_arr, gain_arr, ubs_arr):
            ubs_um = []
            for i in ubs:
                ubs_um.append(lut[int(i)])

            # Assign the first bin lower boundary if the ignore_b1 flag is False.
            if ignore_b1 is False:
                if gain == "Droplet":
                    b1_lb = 1.0                             # For low gain (Droplet)
                elif gain == "Aerosol":
                    b1_lb = 0.3                             # For high gain (Aerosol)
                else:
                    raise ValueError("ERROR: Invalid UCASS gain mode")
                ubs_um.insert(0, b1_lb)                     # Insert value at start of list
            elif ignore_b1 is True:
                pass
            else:
                raise ValueError("ERROR: \'ignore_b1\' is not boolean")

            ubs_um = [float(i) for i in ubs_um]  # Convert to floats for analysis
            ubs_um_list.append(ubs_um)

    else:
        raise TypeError("ERROR: \'sua_data\' is of unrecognised type (type is: %s)" % str(type(sua_data)))

    # Compute the actual bin centres now the values are in intelligible units.
    bin_centres_list = []
    for ubs_um in ubs_um_list:
        bin_centres = []
        if centre_type == "Geometric":
            for i in range(len(ubs_um) - 1):                        # Loop through bins
                centre = float(np.sqrt([ubs_um[i]*ubs_um[i+1]]))    # Equation for geometric mean
                bin_centres.append(centre)
        elif centre_type == "Arithmetic":
            for i in range(len(ubs_um) - 1):                        # Loop through bins
                centre = float((ubs_um[i]+ubs_um[i+1])/2)           # Equation for arithmetic mean
                bin_centres.append(centre)
        else:
            raise ValueError("ERROR: Unrecognised mean type")
        bin_centres_list.append(bin_centres)

    if "CYISUAData" in str(type(sua_data)):
        sua_data.bin_centres_dp_um1 = bin_centres_list[0]
        sua_data.bin_centres_dp_um2 = bin_centres_list[1]
        sua_data.bin_bounds_dp_um1 = ubs_um_list[0]
        sua_data.bin_bounds_dp_um2 = ubs_um_list[1]
    else:
        sua_data.bin_centres_dp_um = bin_centres_list[0]             # Final assignment

    return


def sample_volume(sua_data, altitude_type="GPS", sample_area_m2=0.5e-6, airspeed_type="normal"):
    """
    This script will derive the effective sample volume for a particle counter. This method will differ depending on
    which instrument is used to take the data.
    :param sua_data: The input data object for the instrument
    :param altitude_type: Only used for "SUAData". Either GPS altitude or pressure altitude is used
    :param sample_area_m2: Only used for "SUAData". The area of the laser which is counted as the sample are in m^2.
    :return: Nothing, all assignments effect input object.
    """

    # Ensuring there are no problems with the SUA data class and importing.
    if sua_data.sample_volume_m3 is not None:
        warnings.warn("WARNING: Overwriting existing Analysis")

    # For UCASS or SuperSonde, the sample volume is computed from the sample area and the distance that sample volume
    # has travelled in a time-step.
    if ("SUAData" in str(type(sua_data))) and ("CYI" not in str(type(sua_data))) and ("FMI" not in str(type(sua_data))):

        # Importing variables into namespace
        try:
            gps_alt = sua_data.alt
            gps_alt = gps_alt.astype(float)     # Convert np.array to float types for analysis
        except NameError:
            raise NameError("ERROR: Problem with SUA data object")

        # Computing altitude
        if altitude_type == "Pressure":             # Pressure will become implemented when Temp is recalibrated
            raise NotImplementedError("ERROR: Function not yet supported")
        elif altitude_type == "GPS":
            alt = np.true_divide(gps_alt, 1000)     # Convert GPS altitude to m from mm (default from Pixhawk)
        else:
            raise ValueError("ERROR: Unrecognised altitude analysis type")

        # Compute sample volume
        alt = alt - alt[0]                                      # Normalize altitude to ground level
        integration_length = abs(np.diff(alt, axis=0))          # differentiate to get velocity
        sample_volume_m3 = integration_length*sample_area_m2    # Times by sample area to get volume
        sample_volume_m3 = np.vstack((1, sample_volume_m3))     # Alter length so variable is accepted into class

    # For the static CAS data, the sample volume can be found from total number concentration and the raw counts per bin
    elif "StaticCASData" in str(type(sua_data)):
        try:
            num_conc = sua_data.number_concentration
            counts = sua_data.raw_counts
        except NameError:
            raise NameError("ERROR: Problem with SUA data object")

        c_sum = np.sum(counts, axis=1)          # Get cumulative sum of counts
        size = c_sum.shape                      # Length of loop
        sample_volume_m3 = np.zeros(size)       # Pre allocate array
        index = 0
        for i in c_sum:
            if num_conc[index] == 0:            # Don't divide by 0
                sample_volume_m3[index] = 0
            else:
                sample_volume_m3[index] = i / (num_conc[index] * 1000000)
            index += 1

    elif "StaticFSSPData" in str(type(sua_data)):
        try:
            airspeed = sua_data.airspeed
            sample_area_m2 = sua_data.sample_area_mm2*10**(-6)
            time = sua_data.time
        except NameError:
            raise NameError("ERROR: Problem with SUA data object")

        integration_length = np.diff(time, axis=0)
        integration_length = np.vstack((1, integration_length))
        sample_distance = np.multiply(integration_length, airspeed)
        sample_volume_m3 = np.multiply(sample_distance, sample_area_m2)

    elif ("CYISUAData" in str(type(sua_data))) or ("FMISUAData" in str(type(sua_data))):

        # Importing variables into namespace
        try:
            if airspeed_type == "normal":
                airspeed = sua_data.vz_cms
            elif airspeed_type == "adjusted":
                airspeed = sua_data.adjusted_airspeed
            else:
                raise ValueError("ERROR: Invalid airspeed type")
            airspeed = airspeed.astype(float)  # Convert np.array to float types for analysis
            if "CYI" in str(type(sua_data)):
                period = sua_data.opc_aux1[:, 0]
            else:
                period = sua_data.opc_aux[:, 0]
        except NameError:
            raise NameError("ERROR: Problem with SUA data object")

        airspeed = np.true_divide(airspeed, 100)

        # Compute sample volume
        integration_time = np.multiply(1 / (32.768 * 1000.0), period)
        integration_time = integration_time[None].T
        integration_length = np.multiply(integration_time, airspeed)  # differentiate to get velocity
        sample_volume_m3 = np.multiply(integration_length, sample_area_m2)  # Times by sample area to get volume

    else:
        raise TypeError("ERROR: \'sua_data\' is of unrecognised type (type is: %s)" % str(type(sua_data)))

    sua_data.sample_volume_m3 = sample_volume_m3

    return


def mass_concentration_kgm3(sua_data, material="Water"):
    """
    This function will compute the mass concentration using sample volume and mean bin diameters.
    :param sua_data: input data object from importer.py
    :param material: The measured material, used for density, water by default
    :return: Nothing, all assignments effect input object.
    """

    # Ensuring there are no problems with the SUA data class and importing.
    if "CYISUAData" in str(type(sua_data)):
        if (sua_data.bin_centres_dp_um1 is None) or (sua_data.bin_centres_dp_um1 is None):
            print("INFO: Running bin centre computation")
            bin_centre_dp_um(sua_data)
        if sua_data.sample_volume_m3 is None:
            print("INFO: Running sample volume computation")
            sample_volume(sua_data)
        try:
            bin_centres_arr = [sua_data.bin_centres_dp_um1, sua_data.bin_centres_dp_um2]
            sample_volume_m3 = sua_data.sample_volume_m3
            counts_arr = [sua_data.raw_counts1, sua_data.raw_counts2]
        except NameError:
            raise NameError("ERROR: Problem with SUA data object")
    else:
        if sua_data.bin_centres_dp_um is None:
            print("INFO: Running bin centre computation")
            bin_centre_dp_um(sua_data)
        if sua_data.sample_volume_m3 is None:
            print("INFO: Running sample volume computation")
            sample_volume(sua_data)
        try:
            bin_centres_arr = [sua_data.bin_centres_dp_um]
            sample_volume_m3 = sua_data.sample_volume_m3
            counts_arr = [sua_data.raw_counts]
        except NameError:
            raise NameError("ERROR: Problem with SUA data object")

    # Computing density from 'material' variable and the 'density.txt' file
    module_path = ospath.dirname(ospath.realpath(__file__))
    density_path = module_path + "/density.txt"                 # Path to density file
    if ospath.exists(density_path):
        pass
    else:
        raise ValueError("ERROR: Density path does not exist")
    densities = common.file_to_dict(density_path)               # Convert to LUT
    density = densities[material]                               # Find density from key material

    # Loops for computing mass concentration are below
    mass_conc_buf_list = []
    for counts, bin_centres in zip(counts_arr, bin_centres_arr):
        size = counts.shape
        mass_conc_buf = np.zeros([size[0], 1])
        for i in range(size[0]):                                    # Looping thru altitude/time
            p_vol = []
            for j in range(size[1]):                                # Looping thru bins
                vol_buf = 4.0/3.0 * np.pi * ((bin_centres[j]*10.0**(-6.0))/2.0)**3.0 * float(counts[i, j])
                p_vol.append(vol_buf)
            if sample_volume_m3[i, 0] == 0:                         # Don't divide by 0
                mass_conc_buf[i, 0] = 0
            else:
                mass_conc_buf[i, 0] = float(sum(p_vol)) * density / sample_volume_m3[i, 0]

        mass_conc_buf_list.append(mass_conc_buf)

    if "CYISUAData" in str(type(sua_data)):
        sua_data.mass_concentration1 = mass_conc_buf_list[0]
        sua_data.mass_concentration2 = mass_conc_buf_list[1]
    else:
        sua_data.mass_concentration = mass_conc_buf_list[0]

    return


def num_concentration_m3(sua_data):

    # Ensuring there are no problems with the SUA data class and importing.
    counts = None
    counts1 = None
    counts2 = None
    if sua_data.sample_volume_m3 is None:
        print("INFO: Running sample volume computation")
        sample_volume(sua_data)
    try:
        sample_volume_m3 = sua_data.sample_volume_m3
        if "CYISUAData" in str(type(sua_data)):
            counts1 = sua_data.raw_counts1
            counts2 = sua_data.raw_counts2
        else:
            counts = sua_data.raw_counts
    except NameError:
        raise NameError("ERROR: Problem with SUA data object")

    if "CYISUAData" in str(type(sua_data)):
        size1 = counts1.shape
        size2 = counts2.shape
        num_conc_buf1 = np.zeros([size1[0], 1])
        for i in range(size1[0]):
            if sample_volume_m3[i, 0] == 0:
                num_conc_buf1[i, 0] = 0
            else:
                num_conc_buf1[i, 0] = float(sum(counts1[i, :])) / sample_volume_m3[i, 0]
        sua_data.number_concentration1 = num_conc_buf1
        num_conc_buf2 = np.zeros([size2[0], 1])
        for i in range(size2[0]):
            if sample_volume_m3[i, 0] == 0:
                num_conc_buf2[i, 0] = 0
            else:
                num_conc_buf2[i, 0] = float(sum(counts2[i, :])) / sample_volume_m3[i, 0]
        sua_data.number_concentration2 = num_conc_buf2

        return

    else:
        size = counts.shape
        num_conc_buf = np.zeros([size[0], 1])
        for i in range(size[0]):
            if sample_volume_m3[i, 0] == 0:
                num_conc_buf[i, 0] = 0
            else:
                num_conc_buf[i, 0] = float(sum(counts[i, :])) / sample_volume_m3[i, 0]

        sua_data.number_concentration = num_conc_buf

        return


def dn_dlogdp(sua_data):

    # Ensuring there are no problems with the SUA data class and importing.
    bin_bounds = None
    bin_bounds_arr = None
    if sua_data.sample_volume_m3 is None:
        print("INFO: Running sample volume computation")
        sample_volume(sua_data)
    try:
        counts_arr = None
        counts = None
        sample_volume_array = sua_data.sample_volume_m3
        if ("SUAData" in str(type(sua_data))) and ("CYI" not in str(type(sua_data))):
            _keys = sua_data.alt
            counts = sua_data.raw_counts
            bin_bounds = sua_data.bin_bounds_dp_um

        elif "CYISUAData" in str(type(sua_data)):
            _keys = sua_data.alt
            bin_bounds_arr = [sua_data.bin_bounds_dp_um1, sua_data.bin_bounds_dp_um2]
            counts_arr = [sua_data.raw_counts1, sua_data.raw_counts2]

        elif ("StaticCASData" in str(type(sua_data))) or ("StaticFSSPData" in str(type(sua_data))):
            _keys = sua_data.time
            bin_bounds = sua_data.bins

        else:
            raise TypeError("ERROR: \'sua_data\' is of unrecognised type (type is: %s)" % str(type(sua_data)))

    except NameError:
        raise NameError("ERROR: Problem with SUA data object")

    arr_length = _keys.shape[0]
    keys = [(alt[0], i) for alt, i in zip(_keys, range(arr_length))]

    if "CYISUAData" in str(type(sua_data)):
        dn_dlogdp_dict_arr = []
        for bin_bounds, counts in zip(bin_bounds_arr, counts_arr):
            bins = counts.shape[1]
            dn_dlogdp_dict = {}
            for i in range(arr_length):
                sample_volume_cm3 = sample_volume_array[i] * 1000000.0
                counts_at_key = counts[i, :]
                key = keys[i]
                dn_dlogdp_at_key = []
                for j in range(bins):
                    if sample_volume_cm3 == 0:
                        dn = 0
                    else:
                        dn = counts_at_key[j] / sample_volume_cm3
                    dpl = float(bin_bounds[j]) / 10000.0
                    dpu = float(bin_bounds[j + 1]) / 10000.0
                    dn_dlogdp_in_bin = dn / (np.log10(dpu) - np.log10(dpl))
                    if np.isscalar(dn_dlogdp_in_bin):
                        dn_dlogdp_at_key.append(dn_dlogdp_in_bin)
                    else:
                        dn_dlogdp_at_key.append(dn_dlogdp_in_bin[0])
                    dn_dlogdp_dict[key] = dn_dlogdp_at_key

            dn_dlogdp_dict_arr.append(dn_dlogdp_dict)

        sua_data.dn_dlogdp1 = dn_dlogdp_dict_arr[0]
        sua_data.dn_dlogdp2 = dn_dlogdp_dict_arr[1]

        return

    else:
        bins = counts.shape[1]
        dn_dlogdp_dict = {}
        for i in range(arr_length):
            sample_volume_cm3 = sample_volume_array[i] * 1000000.0
            counts_at_key = counts[i, :]
            key = keys[i]
            dn_dlogdp_at_key = []
            for j in range(bins):
                if sample_volume_cm3 == 0:
                    dn = 0
                else:
                    dn = counts_at_key[j] / sample_volume_cm3
                dpl = float(bin_bounds[j]) / 10000.0
                dpu = float(bin_bounds[j+1]) / 10000.0
                dn_dlogdp_in_bin = dn / (np.log10(dpu) - np.log10(dpl))
                if np.isscalar(dn_dlogdp_in_bin):
                    dn_dlogdp_at_key.append(dn_dlogdp_in_bin)
                else:
                    dn_dlogdp_at_key.append(dn_dlogdp_in_bin[0])
                dn_dlogdp_dict[key] = dn_dlogdp_at_key

        sua_data.dn_dlogdp = dn_dlogdp_dict

        return
