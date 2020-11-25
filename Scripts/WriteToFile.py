from AirborneParticleAnalysis import level1to2
import numpy as np


if __name__ == "__main__":
    data_file = "C:\\Users\\jg17acv\\University of Hertfordshire\\AeroSAM - Documents\\Data\\2020\\20-03-10\\level_1" \
                "\\CYI-FW-UCASS-X2_20200310_115555_00.pdat"
    out_file = "C:\\Users\\jg17acv\\University of Hertfordshire\\AeroSAM - Documents\\Data\\2020\\20-03-10\\level_2" \
               "\\115555_up_nconc_cm3.txt"
    data = level1to2.import_level1(data_file)

    mask = data.up_profile_mask
    alt = data.alt[np.where(mask[:, 0] == 1)]
    n_conc = np.true_divide(data.number_concentration1, 1e6)[np.where(mask[:, 0] == 1)]
    with open(out_file, "w+") as f:
        f.writelines("alt (m),n_conc (#/cm3)\n")
        for a, c in zip(alt, n_conc):
            line_str = str(a[0]) + ',' + str(c[0]) + '\n'
            f.writelines(line_str)


