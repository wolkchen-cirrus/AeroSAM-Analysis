from AirborneParticleAnalysis import common
import os
import datetime
import numpy as np
import warnings


class StaticFSSPData(object):

    def __init__(self, level0_path=None):

        print("INFO: Importing Static FSSP Data")

        # Protected variables to store property data for auxiliary (non-columnated) data.
        self._num_lines = None          # Number of lines
        self._path = None               # Data path
        self._bins = None               # Bin boundaries (upper)
        self._epoch = None              # GPS epoch time
        self._datetime = None           # Human time
        self._row = None                # Row data (used in loop, not for analysis)
        self._row_index = 0             # Row index (ditto)
        self._alt = None                # Altitude ASL in cm
        self._tags = None
        self._level_indicator = None
        self._sample_area_mm2 = None

        # Protected variables to store property data for columnated data (level 0)
        self._time = None                   # Time (epoch) of the line
        self._raw_counts = None             # Raw OPC binned particle counts
        self._number_concentration = None
        self._airspeed = None

        # Protected variables to store property data after level 1 analysis
        self._mass_concentration = None
        self._bin_centres_dp_um = None
        self._dn_dlogdp = None
        self._sample_volume_m3 = None

        if not level0_path:                                     # reading path string
            self.path = common.read_setting("FSSP_level0_data_path")
        else:
            self.path = level0_path

        self.num_lines = common.line_nums(self.path)  # getting number of lines in .csv
        with open(self.path) as f:  # Opening file

            # Assigning auxiliary data to properties.
            lines = f.readlines()
            b0_lb = float(lines[28].split(',')[0].split(">")[1])
            bin_ubs = lines[29].split(',')[1:31]
            bin_ubs = [float(i) for i in bin_ubs]
            bin_ubs.insert(0, b0_lb)
            self.bins = bin_ubs

            self.sample_area_mm2 = 0.414

            self.tags = self.path.split("\\")[-1]
            fssp_date = self.path.split("\\")[-1].split("_")[-2].split(".")[0]
            fssp_y = fssp_date[0:4]
            fssp_m = fssp_date[4:6]
            fssp_d = fssp_date[6:8]
            cas_datetime = \
                datetime.datetime.strptime(fssp_y + "/" + fssp_m + "/" + fssp_d + " 00:00:00", "%Y/%m/%d %H:%M:%S")
            self.epoch = common.utc_to_epoch(cas_datetime)
            self.datetime = self.path.split("\\")[-1].split("_")[-2]

            # Assigning columnated data to properties in loop.
            for i in lines:                         # Loop through lines
                try:
                    self.row = i.split(',')         # Perform row property check
                except (ValueError, TypeError):     # Raised if row is a header/AUX
                    print "INFO: Skipping Row"
                    continue                        # Skip the iteration

                print "INFO: Processing row number %s" % self.row_index

                self.time = float(self.row[0])
                self.airspeed = float(self.row[30])
                self.raw_counts = self.row[31:61]
                self.number_concentration = float(self.row[26])

            self.num_lines = self.row_index

    # These are descriptor objects following the format described in common. The format is general so all the
    # column data is stored under the same conditions, without polluting the namespace of the class.
    time = common.ColumnProperty("time")                    # Time Data
    raw_counts = common.ColumnProperty("raw_counts")        # Raw OPC counts for bins 0-15
    number_concentration = common.ColumnProperty("number_concentration")
    airspeed = common.ColumnProperty("airspeed")

    # These are similar to above but added after the initial import.
    mass_concentration = common.AddedColumn("mass_concentration")
    sample_volume_m3 = common.AddedColumn("sample_volume_m3")

    def check_level(self):
        level_bool = []
        if self.sample_volume_m3 is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.dn_dlogdp is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.bin_centres_dp_um is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)

        self.level_indicator = min(level_bool)

    # The properties that follow are designed to stop the mis-assignment of the AUX values with the data:
    @property
    def sample_area_mm2(self):
        return self._sample_area_mm2

    @sample_area_mm2.setter
    def sample_area_mm2(self, value):
        if not isinstance(value, float):
            raise TypeError
        self._sample_area_mm2 = value

    @property
    def level_indicator(self):
        return self._level_indicator

    @level_indicator.setter
    def level_indicator(self, value):
        if not isinstance(value, int):
            raise TypeError
        self._level_indicator = value

    @property
    def dn_dlogdp(self):
        return self._dn_dlogdp

    @dn_dlogdp.setter
    def dn_dlogdp(self, value):
        if not isinstance(value, dict):
            raise TypeError
        self._dn_dlogdp = value

    @property
    def bin_centres_dp_um(self):
        return self._bin_centres_dp_um

    @bin_centres_dp_um.setter
    def bin_centres_dp_um(self, value):
        if not isinstance(value, list):
            raise TypeError
        self._bin_centres_dp_um = value

    @property
    def num_lines(self):
        return self._num_lines

    @num_lines.setter
    def num_lines(self, value):
        if not isinstance(value, int):
            raise TypeError
        self._num_lines = value

    @property
    def bins(self):
        return self._bins

    @bins.setter
    def bins(self, value):
        if not isinstance(value, list):
            raise TypeError
        if len(value) != 31:
            raise ValueError("ERROR: Must be 31 bin boundaries for FSSP")
        self._bins = value

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, value):
        if not isinstance(value, str):
            raise TypeError
        if not os.path.exists(value) and "level_0" in value:
            raise ValueError("ERROR: Path does not exist")
        if "2018" in value:
            raise ValueError("ERROR: Script only valid for data after 2019")
        self._path = value

    @property
    def epoch(self):
        return self._epoch

    @epoch.setter
    def epoch(self, value):
        if not isinstance(value, int):
            try:
                value = float(value)
            except TypeError:
                raise TypeError("ERROR: Invalid type for epoch")
        self._epoch = value

    @property
    def datetime(self):
        return self._datetime

    @datetime.setter
    def datetime(self, value):
        if isinstance(value, str):
            self._datetime = value
        else:
            raise TypeError("ERROR: Invalid type for datetime")

    @property
    def row(self):
        return self._row

    @row.setter
    def row(self, value):
        value = filter(None, value)
        if not isinstance(value, list):
            raise TypeError("ERROR: Invalid Type for row")
        if len(value) is not 65:
            raise ValueError("ERROR: Must be 65 column rows")
        if isinstance(value[0], str):
            try:
                float(value[0])
            except TypeError:
                raise TypeError("ERROR: Invalid type within list for row")
        self.row_index += 1
        self._row = value

    @property
    def row_index(self):
        return self._row_index

    @row_index.setter
    def row_index(self, value):
        if not isinstance(value, int):
            raise TypeError
        self._row_index = value

    @property
    def tags(self):
        return self._tags

    @tags.setter
    def tags(self, value):
        if value is None:
            print "INFO: No tags assigned"
            return
        elif not isinstance(value, str):
            raise TypeError("ERROR: Tags for CAS must be strings delimited with _")
        else:
            tag_arr = value.split("_")
            module_path = os.path.dirname(os.path.realpath(__file__))
            valid_tags_path = module_path + "/valid_tags.txt"
            with open(valid_tags_path) as f:
                valid_tags = f.read().split(',')
            accepted_tags = []
            for tag in tag_arr:
                try:
                    float(tag.split('.')[0])
                except ValueError:
                    accepted_tags.append(tag)
            for tag in accepted_tags:
                if tag not in valid_tags:
                    warnings.warn("WARNING: Tag %s not in valid tags, check spelling" % tag)
            self._tags = accepted_tags


class StaticCASData(object):

    def __init__(self, level0_path=None):

        print("INFO: Importing Static CAS Data")

        # Protected variables to store property data for auxiliary (non-columnated) data.
        self._num_lines = None          # Number of lines
        self._path = None               # Data path
        self._bins = None               # Bin boundaries (upper)
        self._epoch = None              # GPS epoch time
        self._datetime = None           # Human time
        self._row = None                # Row data (used in loop, not for analysis)
        self._row_index = 0             # Row index (ditto)
        self._alt = None                # Altitude ASL in cm
        self._tags = None
        self._level_indicator = None

        # Protected variables to store property data for columnated data (level 0)
        self._time = None                   # Time (epoch) of the line
        self._raw_counts = None             # Raw OPC binned particle counts
        self._number_concentration = None
        self._lwc_cas_gcm3 = None

        # Protected variables to store property data after level 1 analysis
        self._mass_concentration = None
        self._bin_centres_dp_um = None
        self._dn_dlogdp = None
        self._sample_volume_m3 = None

        if not level0_path:                                     # reading path string
            self.path = common.read_setting("CAS_level0_data_path")
        else:
            self.path = level0_path

        self.num_lines = common.line_nums(self.path)  # getting number of lines in .csv
        with open(self.path) as f:  # Opening file

            # Assigning auxiliary data to properties.
            lines = f.readlines()
            b0_lb = float(lines[24].split("=")[1])
            b0_ub = float(lines[28].split(',')[0].split('>')[1])
            bin_ubs = lines[28].split(',')[1:]
            bin_ubs = [float(i) for i in bin_ubs]
            bin_ubs.insert(0, b0_ub)
            bin_ubs.insert(0, b0_lb)
            self.bins = bin_ubs

            self.tags = self.path.split("\\")[-1]
            cas_date = self.path.split("\\")[-1].split("_")[-2].split(".")[0]
            cas_y = cas_date[0:4]
            cas_m = cas_date[4:6]
            cas_d = cas_date[6:8]
            cas_datetime = \
                datetime.datetime.strptime(cas_y + "/" + cas_m + "/" + cas_d + " 00:00:00", "%Y/%m/%d %H:%M:%S")
            self.epoch = common.utc_to_epoch(cas_datetime)
            self.datetime = self.path.split("\\")[-1].split("_")[-2]

            # Assigning columnated data to properties in loop.
            for i in lines:                         # Loop through lines
                try:
                    self.row = i.split(',')         # Perform row property check
                except (ValueError, TypeError):     # Raised if row is a header/AUX
                    print "INFO: Skipping Row"
                    continue                        # Skip the iteration

                print "INFO: Processing row number %s" % self.row_index

                self.time = float(self.row[0])
                self.raw_counts = self.row[63:93]
                self.number_concentration = float(self.row[53])
                self.lwc_cas_gcm3 = float(self.row[54])

            self.num_lines = self.row_index

    # These are descriptor objects following the format described in common. The format is general so all the
    # column data is stored under the same conditions, without polluting the namespace of the class.
    time = common.ColumnProperty("time")                    # Time Data
    raw_counts = common.ColumnProperty("raw_counts")        # Raw OPC counts for bins 0-15
    number_concentration = common.ColumnProperty("number_concentration")
    lwc_cas_gcm3 = common.ColumnProperty("lwc_cas_gcm3")

    # These are similar to above but added after the initial import.
    mass_concentration = common.AddedColumn("mass_concentration")
    sample_volume_m3 = common.AddedColumn("sample_volume_m3")

    def check_level(self):
        level_bool = []
        if self.sample_volume_m3 is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.dn_dlogdp is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.bin_centres_dp_um is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)

        self.level_indicator = min(level_bool)

    # The properties that follow are designed to stop the mis-assignment of the AUX values with the data:
    @property
    def level_indicator(self):
        return self._level_indicator

    @level_indicator.setter
    def level_indicator(self, value):
        if not isinstance(value, int):
            raise TypeError
        self._level_indicator = value

    @property
    def dn_dlogdp(self):
        return self._dn_dlogdp

    @dn_dlogdp.setter
    def dn_dlogdp(self, value):
        if not isinstance(value, dict):
            raise TypeError
        self._dn_dlogdp = value

    @property
    def bin_centres_dp_um(self):
        return self._bin_centres_dp_um

    @bin_centres_dp_um.setter
    def bin_centres_dp_um(self, value):
        if not isinstance(value, list):
            raise TypeError
        self._bin_centres_dp_um = value

    @property
    def num_lines(self):
        return self._num_lines

    @num_lines.setter
    def num_lines(self, value):
        if not isinstance(value, int):
            raise TypeError
        self._num_lines = value

    @property
    def bins(self):
        return self._bins

    @bins.setter
    def bins(self, value):
        if not isinstance(value, list):
            raise TypeError
        if len(value) != 31:
            raise ValueError("ERROR: Must be 31 bin boundaries for CAS")
        self._bins = value

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, value):
        if not isinstance(value, str):
            raise TypeError
        if not os.path.exists(value) and "level_0" in value:
            raise ValueError("ERROR: Path does not exist")
        if "2018" in value:
            raise ValueError("ERROR: Script only valid for data after 2019")
        self._path = value

    @property
    def epoch(self):
        return self._epoch

    @epoch.setter
    def epoch(self, value):
        if not isinstance(value, int):
            try:
                value = float(value)
            except TypeError:
                raise TypeError("ERROR: Invalid type for epoch")
        self._epoch = value

    @property
    def datetime(self):
        return self._datetime

    @datetime.setter
    def datetime(self, value):
        if isinstance(value, str):
            self._datetime = value
        else:
            raise TypeError("ERROR: Invalid type for datetime")

    @property
    def row(self):
        return self._row

    @row.setter
    def row(self, value):
        value = filter(None, value)
        if not isinstance(value, list):
            raise TypeError("ERROR: Invalid Type for row")
        if len(value) is not 124:
            raise ValueError("ERROR: Must be 98 column rows")
        if isinstance(value[0], str):
            try:
                float(value[0])
            except TypeError:
                raise TypeError("ERROR: Invalid type within list for row")
        self.row_index += 1
        self._row = value

    @property
    def row_index(self):
        return self._row_index

    @row_index.setter
    def row_index(self, value):
        if not isinstance(value, int):
            raise TypeError
        self._row_index = value

    @property
    def tags(self):
        return self._tags

    @tags.setter
    def tags(self, value):
        if value is None:
            print "INFO: No tags assigned"
            return
        elif not isinstance(value, str):
            raise TypeError("ERROR: Tags for CAS must be strings delimited with _")
        else:
            tag_arr = value.split("_")
            module_path = os.path.dirname(os.path.realpath(__file__))
            valid_tags_path = module_path + "/valid_tags.txt"
            with open(valid_tags_path) as f:
                valid_tags = f.read().split(',')
            accepted_tags = []
            for tag in tag_arr:
                try:
                    float(tag.split('.')[0])
                except ValueError:
                    accepted_tags.append(tag)
            for tag in accepted_tags:
                if tag not in valid_tags:
                    warnings.warn("WARNING: Tag %s not in valid tags, check spelling" % tag)
            self._tags = accepted_tags


class SUAData(object):
    """
    This class stores all the data from a UH-AeroSAM .csv data file into a series of protected storage
    variables using properties to specify what gets stored. Level 0 data is filled out immediately upon instantiation of
    the class, Level 1 and above data is filled out by other functions, and left blank at this stage.
    :param level0_path: Path for the raw data, if left as None the importer will infer the data path from "settings.txt"
    """

    def __init__(self, level0_path=None):

        print("INFO: Importing AeroSAM Data")

        # Protected variables to store property data for auxiliary (non-columnated) data.
        self._num_lines = None          # Number of lines
        self._path = None               # Data path
        self._info_string = None        # UCASS info string
        self._bins = None               # Bin boundaries (upper)
        self._gsc = None                # Gain Scaling Coefficient
        self._id = None                 # UCASS ID
        self._epoch = None              # GPS epoch time
        self._datetime = None           # Human time
        self._trash = None              # Trash flag
        self._tags = None               # User assigned tags for data, used for searching
        self._row = None                # Row data (used in loop, not for analysis)
        self._row_index = 0             # Row index (ditto)
        self._level_indicator = 0

        # Protected variables to store property data for columnated data (level 0)
        self._time = None               # Time (epoch) of the line
        self._press_hpa = None          # Atmospheric pressure in hPa
        self._lat = None                # Latitude co-ordinate
        self._lon = None                # Longitude co-ordinate
        self._alt = None                # Altitude ASL in cm
        self._vz_cms = None             # 'z' velocity in cm/s
        self._temp_deg_c = None         # Temperature in degrees C
        self._rh_true = None            # True (temp-corrected) relative humidity as a %
        self._raw_counts = None         # Raw OPC binned particle counts
        self._m_tof = None              # Mean time of flight data
        self._opc_aux = None            # Auxiliary OPC data (column specific e.g. glitch trap)

        # Protected variables to store property data after level 1 analysis
        self._down_profile_mask = None      # multi-dimensional np.array which masks out waiting time (Downwards)
        self._up_profile_mask = None        # multi-dimensional np.array which masks out waiting time (Upwards)
        self._profile_number = None         # Number of profiles in one data object
        self._ucass_lut_aerosol = None
        self._ucass_lut_droplet = None
        self._mass_concentration = None
        self._number_concentration = None
        self._bin_centres_dp_um = None
        self._sample_volume_m3 = None
        self._bin_bounds_dp_um = None
        self._dn_dlogdp = None

        # Recording the file data to class properties. The data path is specified in the settings.txt file. This will
        # start by getting the AUX data, then move onto the columnated data in a loop.
        if not level0_path:                                     # reading path string
            self.path = common.read_setting("SAM_level0_data_path")
        else:
            self.path = level0_path
        self.num_lines = common.line_nums(self.path)            # getting number of lines in .csv
        filename_time = self.path.split("\\")[-1].split("_")[-2]
        filename_date = self.path.split("\\")[-1].split("_")[-3]
        with open(self.path) as f:                              # Opening file
            lines = f.readlines()

            # Assigning auxiliary data to properties.
            self.info_string = lines[1].split(',')[0]           # Getting info string
            self.bins = lines[3].split(',')[:16]                # Getting bin ADC boundaries
            self.gsc = lines[3].split(',')[16]                  # Getting gain scaling coefficient
            if lines[3].split(',')[17] is ")" or "(":           # Accounting for errors in older data
                self.id = lines[3].split(',')[18]               # If yes skip it
            else:
                self.id = lines[3].split(',')[17]               # If no read the next one
            self.epoch = lines[0].split(',')[2]                 # Getting GPS epoch time
            self.datetime = filename_date + filename_time       # Converting to human date
            self.trash = lines[0].split(',')[3]                 # Getting trash indicator
            self.tags = lines[0].split(',')[4]                  # Assigning tags for data

            # Assigning columnated data to properties in loop.
            for i in lines:                                     # Loop through lines
                try:
                    self.row = i.split(',')                     # Perform row property check
                except (ValueError, TypeError):                 # Raised if row is a header/AUX
                    print "INFO: Skipping Row"
                    continue                                    # Skip the iteration

                # Divide up the row property attribute, and append to the column properties. Note that the appending is
                # done automatically with in common.ColumnProperty() class when the __set__ method is called upon the
                # assignment of an attribute.
                self.time = self.row[0]
                self.press_hpa = self.row[1]
                self.lat = self.row[2]
                self.lon = self.row[3]
                self.alt = self.row[4]
                self.vz_cms = self.row[5]
                self.temp_deg_c = self.row[6]
                self.rh_true = self.row[7]
                self.raw_counts = self.row[8:24]
                self.m_tof = self.row[24:28]
                self.opc_aux = self.row[28:]

            self.num_lines = self.row_index

    # These are descriptor objects following the format described in common. The format is general so all the
    # column data is stored under the same conditions, without polluting the namespace of the class.
    time = common.ColumnProperty("time")                    # Time Data
    press_hpa = common.ColumnProperty("press_hpa")          # Pressure in hPa
    lat = common.ColumnProperty("lat")                      # Latitude co-ordinate
    lon = common.ColumnProperty("lon")                      # Longitude co-ordinate
    alt = common.ColumnProperty("alt")                      # Altitude in mm
    vz_cms = common.ColumnProperty("vz_cms")                # Vertical velocity in cm/s
    temp_deg_c = common.ColumnProperty("temp_deg_c")        # Temperature in degrees C
    rh_true = common.ColumnProperty("rh_true")              # True relative humidity as %
    raw_counts = common.ColumnProperty("raw_counts")        # Raw OPC counts for bins 0-15
    m_tof = common.ColumnProperty("m_tof")                  # Mean time of flight data
    opc_aux = common.ColumnProperty("opc_aux")              # Auxiliary OPC data (glitch trap etc.)

    # These are similar to above but added after the initial import.
    sample_volume_m3 = common.AddedColumn("sample_volume_m3")
    mass_concentration = common.AddedColumn("mass_concentration")
    number_concentration = common.AddedColumn("number_concentration")

    def check_level(self):
        level_bool = []
        if self.sample_volume_m3 is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.dn_dlogdp is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.bin_centres_dp_um is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.bin_bounds_dp_um is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.mass_concentration is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.number_concentration is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.up_profile_mask is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if (self.ucass_lut_aerosol is not None) or (self.ucass_lut_droplet is not None):
            level_bool.append(1)
        else:
            level_bool.append(0)

        self.level_indicator = min(level_bool)

    # The properties that follow are designed to stop the mis-assignment of the AUX values with the data:
    @property
    def level_indicator(self):
        return self._level_indicator

    @level_indicator.setter
    def level_indicator(self, value):
        if not isinstance(value, int):
            raise TypeError
        self._level_indicator = value

    @property
    def dn_dlogdp(self):
        return self._dn_dlogdp

    @dn_dlogdp.setter
    def dn_dlogdp(self, value):
        if not isinstance(value, dict):
            raise TypeError
        self._dn_dlogdp = value

    @property
    def bin_bounds_dp_um(self):
        return self._bin_bounds_dp_um

    @bin_bounds_dp_um.setter
    def bin_bounds_dp_um(self, value):
        if not isinstance(value, list):
            raise TypeError
        self._bin_bounds_dp_um = value

    @property
    def bin_centres_dp_um(self):
        return self._bin_centres_dp_um

    @bin_centres_dp_um.setter
    def bin_centres_dp_um(self, value):
        if not isinstance(value, list):
            raise TypeError
        self._bin_centres_dp_um = value

    @property
    def ucass_lut_aerosol(self):
        return self._ucass_lut_aerosol

    @ucass_lut_aerosol.setter
    def ucass_lut_aerosol(self, value):
        if isinstance(value, dict):
            self._ucass_lut_aerosol = value
        else:
            raise TypeError

    @property
    def ucass_lut_droplet(self):
        return self._ucass_lut_droplet

    @ucass_lut_droplet.setter
    def ucass_lut_droplet(self, value):
        if isinstance(value, dict):
            self._ucass_lut_droplet = value
        else:
            raise TypeError

    @property
    def tags(self):
        return self._tags

    @tags.setter
    def tags(self, value):
        if value is None:
            print "INFO: No tags assigned"
            return
        elif not isinstance(value, str):
            raise TypeError("ERROR: Tags must be strings delimited with |")
        else:
            tag_arr = value.split("|")
            module_path = os.path.dirname(os.path.realpath(__file__))
            valid_tags_path = module_path + "/valid_tags.txt"
            with open(valid_tags_path) as f:
                valid_tags = f.read().split(',')
            for tag in tag_arr:
                if tag not in valid_tags:
                    warnings.warn("WARNING: Tag %s not in valid tags, check spelling" % tag)
            self._tags = tag_arr

    @property
    def up_profile_mask(self):
        return self._up_profile_mask

    @up_profile_mask.setter
    def up_profile_mask(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError("ERROR: Profile mask must be ndarray")
        self.profile_number = int(value.shape[1])
        self._up_profile_mask = value

    @property
    def down_profile_mask(self):
        return self._down_profile_mask

    @down_profile_mask.setter
    def down_profile_mask(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError("ERROR: Profile mask must be ndarray")
        self.profile_number = int(value.shape[1])
        self._down_profile_mask = value

    @property
    def profile_number(self):
        return self._profile_number

    @profile_number.setter
    def profile_number(self, value):
        if not isinstance(value, int):
            raise TypeError
        self._profile_number = value * 2

    @property
    def num_lines(self):
        return self._num_lines

    @num_lines.setter
    def num_lines(self, value):
        if not isinstance(value, int):
            raise TypeError
        self._num_lines = value

    @property
    def info_string(self):
        return self._info_string

    @info_string.setter
    def info_string(self, value):
        if "FirmwareVer=" not in value:
            warnings.warn("WARNING: UCASS not connected for flight")
        if not isinstance(value, str):
            raise TypeError
        self._info_string = value

    @property
    def bins(self):
        return self._bins

    @bins.setter
    def bins(self, value):
        if not isinstance(value, list):
            raise TypeError
        if len(value) > 16:
            raise ValueError
        self._bins = value

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, value):
        if not isinstance(value, str):
            raise TypeError
        if not os.path.exists(value) and "level_0" in value:
            raise ValueError("ERROR: Path does not exist")
        if "2018" in value:
            raise ValueError("ERROR: Script only valid for data after 2019")
        self._path = value

    @property
    def gsc(self):
        return self._gsc

    @gsc.setter
    def gsc(self, value):
        if not isinstance(value, float):
            if isinstance(value, str):
                print "INFO: Converting GSC to float"
                if '(' or ')' in value:
                    value = value.replace('(', '').replace(')', '')
                try:
                    value = float(value)
                except(TypeError, ValueError):
                    raise TypeError("ERROR: Invalid Type")
            else:
                raise TypeError
        if int(value) is not 1:
            raise ValueError("ERROR: GSC is not equal to 1, double check calibration")
        self._gsc = value

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, value):
        if not isinstance(value, int):
            try:
                value = int(value)
            except TypeError:
                print "ERROR: ID must be integer"
                raise TypeError
        if not 0 <= value <= 255:
            raise ValueError("ERROR: ID must be a 8-bit integer")
        self._id = value

    @property
    def epoch(self):
        return self._epoch

    @epoch.setter
    def epoch(self, value):
        if not isinstance(value, int):
            try:
                value = float(value)
            except TypeError:
                raise TypeError("ERROR: Invalid type for epoch")
        self._epoch = value

    @property
    def datetime(self):
        return self._datetime

    @datetime.setter
    def datetime(self, value):
        if isinstance(value, str):
            self._datetime = value
        else:
            raise TypeError("ERROR: Invalid type for datetime")

    @property
    def trash(self):
        return self._trash

    @trash.setter
    def trash(self, value):
        try:
            value = bool(int(value))
            self._trash = value
            if value is 1:
                warnings.warn("WARNING: Data has been user-specified as trash")
        except TypeError:
            warnings.warn("WARNING: File data \'trash\' boolean not specified, treat with caution.")
            pass

    @property
    def row(self):
        return self._row

    @row.setter
    def row(self, value):
        value = filter(None, value)
        if not isinstance(value, list):
            raise TypeError("ERROR: Invalid Type for row")
        if len(value) is not 33:
            raise ValueError("ERROR: Must be 33 column rows")
        if isinstance(value[0], str):
            try:
                float(value[0])
            except TypeError:
                raise TypeError("ERROR: Invalid type within list for row")
        self.row_index += 1
        self._row = value

    @property
    def row_index(self):
        return self._row_index

    @row_index.setter
    def row_index(self, value):
        if not isinstance(value, int):
            raise TypeError
        self._row_index = value


class CYISUAData(object):
    """
    This class stores all the data from a Cyprus Institute SUA .csv data file into a series of protected storage
    variables using properties to specify what gets stored. Level 0 data is filled out immediately upon instantiation of
    the class, Level 1 and above data is filled out by other functions, and left blank at this stage.
    :param level0_path: Path for the raw data, if left as None the importer will infer the data path from "settings.txt"
    """

    def __init__(self, level0_path=None):

        print("INFO: Importing AeroSAM Data")

        # Protected variables to store property data for auxiliary (non-columnated) data.
        self._num_lines = None          # Number of lines
        self._path = None               # Data path
        self._fd_path = None            # Flight data path
        self._metd_path = None          # MET data path
        self._bins1 = None              # Bin boundaries (upper)
        self._bins2 = None              # Bin boundaries (upper)
        self._epoch = None              # GPS epoch time
        self._datetime = None           # Human time
        self._trash = None              # Trash flag
        self._tags = None               # User assigned tags for data, used for searching
        self._row = None                # Row data (used in loop, not for analysis)
        self._row_index = 0             # Row index (ditto)
        self._level_indicator = 0

        # Protected variables to store property data for columnated data (level 0)
        self._time = None               # Time (epoch) of the line
        self._press_hpa = None          # Atmospheric pressure in hPa
        self._lat = None                # Latitude co-ordinate
        self._lon = None                # Longitude co-ordinate
        self._alt = None                # Altitude ASL in m
        self._vz_cms = None             # 'z' velocity in cm/s
        self._temp_deg_c = None         # Temperature in degrees C
        self._rh_true = None            # True (temp-corrected) relative humidity as a %
        self._raw_counts1 = None        # Raw OPC binned particle counts
        self._m_tof1 = None             # Mean time of flight data
        self._opc_aux1 = None           # Auxiliary OPC data (column specific e.g. glitch trap)
        self._raw_counts2 = None        # Raw OPC binned particle counts
        self._m_tof2 = None             # Mean time of flight data
        self._opc_aux2 = None           # Auxiliary OPC data (column specific e.g. glitch trap)

        # Protected variables to store property data after level 1 analysis
        self._down_profile_mask = None      # multi-dimensional np.array which masks out waiting time (Downwards)
        self._up_profile_mask = None        # multi-dimensional np.array which masks out waiting time (Upwards)
        self._profile_number = None         # Number of profiles in one data object
        self._ucass_lut1 = None
        self._ucass_lut2 = None
        self._ucass_name1 = None
        self._ucass_name2 = None
        self._mass_concentration = None
        self._number_concentration = None
        self._bin_centres_dp_um = None
        self._sample_volume_m3 = None
        self._bin_bounds_dp_um = None
        self._dn_dlogdp = None

        # Recording the file data to class properties. The data path is specified in the settings.txt file. This will
        # start by getting the AUX data, then move onto the columnated data in a loop.
        if not level0_path:                                     # reading path string
            self.path = common.read_setting("CYISUA_level0_data_path")
        else:
            self.path = level0_path
        self.num_lines = common.line_nums(self.path)            # getting number of lines in .csv
        filename_time = self.path.split("\\")[-1].split("_")[-2]
        filename_date = self.path.split("\\")[-1].split("_")[-3]

        test_types = ["CYISUA_FD_path", "CYISUA_METD_path"]
        fd_paths = []
        for fd_type in test_types:
            fd_dir = common.read_setting(fd_type)
            fd_files = os.listdir(fd_dir)
            date_diff = []
            for fd_file in fd_files:
                fd_time = fd_file.split("_")[-2]
                date_diff.append(abs(filename_time-fd_time))
            val = min(date_diff)
            matched_date_index = date_diff.index(val)
            fd_paths.append(fd_dir + "\\" + fd_files[matched_date_index])
        self.fd_path = fd_paths[0]
        self.metd_path = fd_paths[1]

        with open(self.path) as f:                              # Opening file
            lines = f.readlines()

            # Assigning auxiliary data to properties.
            self.bins1 = lines[3].split(',')[:16]               # Getting bin ADC boundaries for UCASS 1
            self.bins2 = lines[3].split(',')[16:32]             # Getting bin ADC boundaries for UCASS 2
            self.datetime = filename_date + filename_time       # Converting to human date
            self.trash = lines[0].split(',')[3]                 # Getting trash indicator
            self.tags = lines[0].split(',')[4]                  # Assigning tags for data

            # Assigning columnated data to properties in loop.
            for i in lines:                                     # Loop through lines
                try:
                    self.row = i.split(',')                     # Perform row property check
                except (ValueError, TypeError):                 # Raised if row is a header/AUX
                    print "INFO: Skipping Row"
                    continue                                    # Skip the iteration

                # Divide up the row property attribute, and append to the column properties. Note that the appending is
                # done automatically with in common.ColumnProperty() class when the __set__ method is called upon the
                # assignment of an attribute.
                self.time = common.hhmmss_to_sec("".join(self.row[6].split(":")))
                self.press_hpa = self.row[1]
                self.lat = self.row[9]
                self.lon = self.row[8]
                self.alt = self.row[10]
                # self.vz_cms = self.row[5]
                # self.temp_deg_c = self.row[6]
                # self.rh_true = self.row[7]
                self.ucass_name1 = self.row[13]
                self.raw_counts1 = self.row[14:30]
                self.m_tof1 = self.row[30:34]
                self.opc_aux1 = self.row[34:39]
                self.ucass_name2 = self.row[39]
                self.raw_counts2 = self.row[40:56]
                self.m_tof2 = self.row[56:60]
                self.opc_aux2 = self.row[60:65]

            self.num_lines = self.row_index

    # These are descriptor objects following the format described in common. The format is general so all the
    # column data is stored under the same conditions, without polluting the namespace of the class.
    time = common.ColumnProperty("time")                    # Time Data
    press_hpa = common.ColumnProperty("press_hpa")          # Pressure in hPa
    lat = common.ColumnProperty("lat")                      # Latitude co-ordinate
    lon = common.ColumnProperty("lon")                      # Longitude co-ordinate
    alt = common.ColumnProperty("alt")                      # Altitude in mm
    vz_cms = common.ColumnProperty("vz_cms")                # Vertical velocity in cm/s
    temp_deg_c = common.ColumnProperty("temp_deg_c")        # Temperature in degrees C
    rh_true = common.ColumnProperty("rh_true")              # True relative humidity as %
    raw_counts1 = common.ColumnProperty("raw_counts1")      # Raw OPC counts for bins 0-15
    m_tof1 = common.ColumnProperty("m_tof1")                # Mean time of flight data
    opc_aux1 = common.ColumnProperty("opc_aux1")            # Auxiliary OPC data (glitch trap etc.)
    raw_counts2 = common.ColumnProperty("raw_counts2")
    m_tof2 = common.ColumnProperty("m_tof2")
    opc_aux2 = common.ColumnProperty("opc_aux2")

    # These are similar to above but added after the initial import.
    sample_volume_m3 = common.AddedColumn("sample_volume_m3")
    mass_concentration = common.AddedColumn("mass_concentration")
    number_concentration = common.AddedColumn("number_concentration")

    def check_level(self):
        level_bool = []
        if self.sample_volume_m3 is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.dn_dlogdp is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.bin_centres_dp_um is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.bin_bounds_dp_um is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.mass_concentration is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.number_concentration is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if self.up_profile_mask is not None:
            level_bool.append(1)
        else:
            level_bool.append(0)
        if (self.ucass_lut1 is not None) or (self.ucass_lut2 is not None):
            level_bool.append(1)
        else:
            level_bool.append(0)

        self.level_indicator = min(level_bool)

    # The properties that follow are designed to stop the mis-assignment of the AUX values with the data:
    @property
    def ucass_name1(self):
        return self._ucass_name1

    @ucass_name1.setter
    def ucass_name1(self, value):
        if not isinstance(value, str):
            raise TypeError
        elif self._ucass_name1 is None:
            self._ucass_name1 = value
            return
        elif self._ucass_name1 != value:
            raise ValueError("ERROR: UCASS ID changes in column")
        else:
            self._ucass_name1 = value
            return

    @property
    def ucass_name2(self):
        return self._ucass_name2

    @ucass_name2.setter
    def ucass_name2(self, value):
        if not isinstance(value, str):
            raise TypeError
        elif self._ucass_name2 is None:
            self._ucass_name2 = value
            return
        elif self._ucass_name2 != value:
            raise ValueError("ERROR: UCASS ID changes in column")
        else:
            self._ucass_name2 = value
            return

    @property
    def level_indicator(self):
        return self._level_indicator

    @level_indicator.setter
    def level_indicator(self, value):
        if not isinstance(value, int):
            raise TypeError
        self._level_indicator = value

    @property
    def dn_dlogdp(self):
        return self._dn_dlogdp

    @dn_dlogdp.setter
    def dn_dlogdp(self, value):
        if not isinstance(value, dict):
            raise TypeError
        self._dn_dlogdp = value

    @property
    def bin_bounds_dp_um(self):
        return self._bin_bounds_dp_um

    @bin_bounds_dp_um.setter
    def bin_bounds_dp_um(self, value):
        if not isinstance(value, list):
            raise TypeError
        self._bin_bounds_dp_um = value

    @property
    def bin_centres_dp_um(self):
        return self._bin_centres_dp_um

    @bin_centres_dp_um.setter
    def bin_centres_dp_um(self, value):
        if not isinstance(value, list):
            raise TypeError
        self._bin_centres_dp_um = value

    @property
    def ucass_lut1(self):
        return self._ucass_lut1

    @ucass_lut1.setter
    def ucass_lut1(self, value):
        if isinstance(value, dict):
            self._ucass_lut1 = value
        else:
            raise TypeError

    @property
    def ucass_lut2(self):
        return self._ucass_lut2

    @ucass_lut2.setter
    def ucass_lut2(self, value):
        if isinstance(value, dict):
            self._ucass_lut2 = value
        else:
            raise TypeError

    @property
    def tags(self):
        return self._tags

    @tags.setter
    def tags(self, value):
        if value is None:
            print "INFO: No tags assigned"
            return
        elif not isinstance(value, str):
            raise TypeError("ERROR: Tags must be strings delimited with |")
        else:
            tag_arr = value.split("|")
            module_path = os.path.dirname(os.path.realpath(__file__))
            valid_tags_path = module_path + "/valid_tags.txt"
            with open(valid_tags_path) as f:
                valid_tags = f.read().split(',')
            for tag in tag_arr:
                if tag not in valid_tags:
                    warnings.warn("WARNING: Tag %s not in valid tags, check spelling" % tag)
            self._tags = tag_arr

    @property
    def up_profile_mask(self):
        return self._up_profile_mask

    @up_profile_mask.setter
    def up_profile_mask(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError("ERROR: Profile mask must be ndarray")
        self.profile_number = int(value.shape[1])
        self._up_profile_mask = value

    @property
    def down_profile_mask(self):
        return self._down_profile_mask

    @down_profile_mask.setter
    def down_profile_mask(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError("ERROR: Profile mask must be ndarray")
        self.profile_number = int(value.shape[1])
        self._down_profile_mask = value

    @property
    def profile_number(self):
        return self._profile_number

    @profile_number.setter
    def profile_number(self, value):
        if not isinstance(value, int):
            raise TypeError
        self._profile_number = value * 2

    @property
    def num_lines(self):
        return self._num_lines

    @num_lines.setter
    def num_lines(self, value):
        if not isinstance(value, int):
            raise TypeError
        self._num_lines = value

    @property
    def bins1(self):
        return self._bins1

    @bins1.setter
    def bins1(self, value):
        if not isinstance(value, list):
            raise TypeError
        if len(value) > 16:
            raise ValueError
        self._bins1 = value

    @property
    def bins2(self):
        return self._bins2

    @bins2.setter
    def bins2(self, value):
        if not isinstance(value, list):
            raise TypeError
        if len(value) > 16:
            raise ValueError
        self._bins2 = value

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, value):
        if not isinstance(value, str):
            raise TypeError
        if not os.path.exists(value) and "level_0" in value:
            raise ValueError("ERROR: Path does not exist")
        if "2018" in value:
            raise ValueError("ERROR: Script only valid for data after 2019")
        self._path = value

    @property
    def fd_path(self):
        return self._fd_path

    @fd_path.setter
    def fd_path(self, value):
        if not isinstance(value, str):
            raise TypeError
        if not os.path.exists(value):
            raise ValueError("ERROR: Path does not exist")
        if "2018" in value:
            raise ValueError("ERROR: Script only valid for data after 2019")
        self._fd_path = value

    @property
    def metd_path(self):
        return self._metd_path

    @metd_path.setter
    def metd_path(self, value):
        if not isinstance(value, str):
            raise TypeError
        if not os.path.exists(value):
            raise ValueError("ERROR: Path does not exist")
        if "2018" in value:
            raise ValueError("ERROR: Script only valid for data after 2019")
        self._metd_path = value

    @property
    def epoch(self):
        return self._epoch

    @epoch.setter
    def epoch(self, value):
        if not isinstance(value, int):
            try:
                value = float(value)
            except TypeError:
                raise TypeError("ERROR: Invalid type for epoch")
        self._epoch = value

    @property
    def datetime(self):
        return self._datetime

    @datetime.setter
    def datetime(self, value):
        if isinstance(value, str):
            self._datetime = value
        else:
            raise TypeError("ERROR: Invalid type for datetime")

    @property
    def trash(self):
        return self._trash

    @trash.setter
    def trash(self, value):
        try:
            value = bool(int(value))
            self._trash = value
            if value is 1:
                warnings.warn("WARNING: Data has been user-specified as trash")
        except TypeError:
            warnings.warn("WARNING: File data \'trash\' boolean not specified, treat with caution.")
            pass

    @property
    def row(self):
        return self._row

    @row.setter
    def row(self, value):
        value = filter(None, value)
        if not isinstance(value, list):
            raise TypeError("ERROR: Invalid Type for row")
        if len(value) is not 65:
            raise ValueError("ERROR: Must be 65 column rows")
        if isinstance(value[0], str):
            try:
                float(value[0])
            except TypeError:
                raise TypeError("ERROR: Invalid type within list for row")
        self.row_index += 1
        self._row = value

    @property
    def row_index(self):
        return self._row_index

    @row_index.setter
    def row_index(self, value):
        if not isinstance(value, int):
            raise TypeError
        self._row_index = value
