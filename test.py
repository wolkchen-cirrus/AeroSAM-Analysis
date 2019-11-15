import importer
import level0to1
import numpy as np


if __name__ == "__main__":
    data = importer.SUAData()
    level0to1.split_by_pressure(data)
    level0to1.assign_ucass_lut(data)
    level0to1.bin_centre_dp_um(data)
    level0to1.ucass_sample_volume(data)
    level0to1.mass_concentration_umm3(data)
    level0to1.num_concentration_m3(data)

    m_conc = data.mass_concentration
    n_conc = data.number_concentration
    split = data.up_profile_mask
    m_conc_prof1 = np.multiply(m_conc, split)
    n_conc_prof1 = np.multiply(n_conc, split)

pass
