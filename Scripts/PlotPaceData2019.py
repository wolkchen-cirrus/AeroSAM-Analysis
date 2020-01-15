from AirborneParticleAnalysis import PacePlots
from AirborneParticleAnalysis import level1to2
from os import listdir


if __name__ == "__main__":
    print "===================== Executing the Plotting Script for the PaCE 2019 Campaign Data =====================\n"

    data_dir = "C:\\Users\\jg17acv\\University of Hertfordshire\\AeroSAM - Documents\\Data\\2019\\19-09-28\\level_1"
    station_altitude_asl_mm = 560000
    data_files = listdir(data_dir)
    data_dict = {}
    for file_name in data_files:
        print "INFO: Importing file - %s" % file_name
        key_name = file_name.split(".")[0]
        file_path = data_dir + "\\" + file_name
        data_dict[key_name] = level1to2.import_level1(file_path)
    dn_dlogdp = []
    times = []

    print "INFO: Analysis pipeline first pass"
    for key in data_dict:
        if ("CAS" in key) or ("FSSP" in key):
            print "INFO: Skipping CAS data and FSSP data on first round"
            continue
        else:
            dn_buf = data_dict[key].dn_dlogdp
            row_number = level1to2.fetch_row(altitude=station_altitude_asl_mm, level1_data=data_dict[key])

