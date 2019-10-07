import common
import os
import datetime
import numpy as np


class SUAData(object):
    """
    This class stores all the data from a UH-AeroSAM .csv data file into a series of protected storage
    variables using properties to specify what gets stored. Level 0 data is filled out immediately upon instantiation of
    the class, Level 1 and above data is filled out by other functions, and left blank at this stage.
    :param level0_path: Path for the raw data, if left as None the importer will infer the data path from "settings.txt"
    """

    def __init__(self, level0_path=None):

        # Protected variables to store property data for auxiliary (non-columnated) data.
        self._num_lines = None          # Number of lines
        self._path = None               # Data path
        self._info_string = None        # UCASS info string
        self._bins = None               # Bin boundaries
        self._gsc = None                # Gain Scaling Coefficient
        self._id = None                 # UCASS ID
        self._epoch = None              # GPS epoch time
        self._datetime = None           # Human time
        self._trash = None              # Trash flag
        self._tags = None               # User assigned tags for data, used for searching
        self._row = None                # Row data (used in loop, not for analysis)
        self._row_index = 0             # Row index (ditto)

        # Protected variables to store property data for columnated data (level 0)
        self._time = None               # Time (epoch) of the line
        self._press_hpa = None          # Atmospheric pressure in hPa
        self._lat = None                # Latitude co-ordinate
        self._lon = None                # Longitude co-ordinate
        self._alt = None                # Altitude ASL in mm
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

        # Recording the file data to class properties. The data path is specified in the settings.txt file. This will
        # start by getting the AUX data, then move onto the columnated data in a loop.
        if not level0_path:                                     # reading path string
            self.path = common.read_setting("level0_data_path")
        else:
            self.path = level0_path
        self.num_lines = common.line_nums(self.path)            # getting number of lines in .csv
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
            self.datetime = self.epoch                          # Converting to human date
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

    # The properties that follow are designed to stop the mis-assignment of the AUX values with the data:
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
            valid_tags_path = common.read_setting("tags_path")
            with open(valid_tags_path) as f:
                valid_tags = f.read().split(',')
            for tag in tag_arr:
                if tag not in valid_tags:
                    print "WARNING: Tag %s not in valid tags, check spelling" % tag
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
            print "WARNING: UCASS not connected for flight"
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
        if not os.path.exists(value):
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
        if isinstance(value, int):
            value = float(value/1000)
            self._datetime = datetime.datetime.fromtimestamp(value).strftime('%Y-%m-%d %H:%M:%S')
        else:
            try:
                value = int(value)
                value = float(value / 1000)
                self._datetime = datetime.datetime.fromtimestamp(value).strftime('%Y-%m-%d %H:%M:%S')
            except TypeError:
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
                print "WARNING: Data has been user-specified as trash"
        except TypeError:
            print "WARNING: File data \'trash\' boolean not specified, treat with caution."
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
