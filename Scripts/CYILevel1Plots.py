from AirborneParticleAnalysis import level1to2
from os import listdir


if __name__ == "__main__":
    print "============ Executing the Plotting Script for CYI SUA Level 1 Data ============\n"

    # Starting the import of the level 1 data, the script "convert_all_0to1.py" must be run first to obtain this level 1
    # data in the first instance.
    data_dir = "C:\\Users\\JGirdwood\\University of Hertfordshire\\AeroSAM - Documents\\Data\\2020\\20-03-10\\level_1"
    data_files = listdir(data_dir)                                                      # Getting files list
    data_dict = {}
    for file_name in data_files:
        print "INFO: Importing file - %s" % file_name
        key_name = file_name.split(".")[0]                                              # Dict key name is file name
        file_path = data_dir + "\\" + file_name                                         # Find full file path
        data_dict[key_name] = level1to2.import_level1(file_path)                        # Importing, uses "cPickle"
