import importer
import level0to1
import numpy as np
from matplotlib import pyplot as ppl


if __name__ == "__main__":

    sam_data = importer.SUAData()
    # cas_data = importer.StaticCASData()

    level0to1.split_by_pressure(sam_data)
    level0to1.assign_ucass_lut(sam_data)
    level0to1.bin_centre_dp_um(sam_data)
    level0to1.sample_volume(sam_data)
    level0to1.mass_concentration_kgm3(sam_data)
    level0to1.num_concentration_m3(sam_data)
    level0to1.dn_dlogdp(sam_data)

    level0to1.export_level1(sam_data)

    # level0to1.bin_centre_dp_um(cas_data)
    # level0to1.sample_volume(cas_data)
    # level0to1.dn_dlogdp(cas_data)

    m_conc = sam_data.mass_concentration
    n_conc = sam_data.number_concentration
    split = sam_data.up_profile_mask
    alt_asl_cm = sam_data.alt
    m_conc_prof1 = np.multiply(m_conc, split)
    n_conc_prof1 = np.multiply(n_conc, split)
    dn_dlogdp = sam_data.dn_dlogdp
    bin_centres_um = sam_data.bin_centres_dp_um
    dn1 = dn_dlogdp[float(alt_asl_cm[60])]

    # ppl.plot(np.multiply(n_conc[55:80]/1000000.0, split[55:80]), np.multiply(alt_asl_cm[55:80], split[55:80]))
    ppl.plot(np.asarray(bin_centres_um), np.asarray(dn_dlogdp[float(alt_asl_cm[70])]))
    ppl.show()

pass
