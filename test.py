import importer
import level0to1
import numpy as np
from matplotlib import pyplot as ppl


if __name__ == "__main__":

    data = importer.SUAData()

    level0to1.split_by_pressure(data)
    level0to1.assign_ucass_lut(data)
    level0to1.bin_centre_dp_um(data)
    level0to1.ucass_sample_volume(data)
    level0to1.mass_concentration_kgm3(data)
    level0to1.num_concentration_m3(data)
    level0to1.dn_dlogdp(data)

    m_conc = data.mass_concentration
    n_conc = data.number_concentration
    split = data.up_profile_mask
    alt_asl_cm = data.alt
    m_conc_prof1 = np.multiply(m_conc, split)
    n_conc_prof1 = np.multiply(n_conc, split)
    dn_dlogdp = data.dn_dlogdp
    bin_centres_um = data.bin_centres_dp_um
    dn1 = dn_dlogdp[float(alt_asl_cm[60])]

    # ppl.plot(np.multiply(n_conc[55:80]/1000000.0, split[55:80]), np.multiply(alt_asl_cm[55:80], split[55:80]))
    ppl.plot(np.asarray(bin_centres_um), np.asarray(dn_dlogdp[float(alt_asl_cm[70])]))
    ppl.show()

pass
