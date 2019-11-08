import importer
import level0to1


data = importer.SUAData()
level0to1.split_by_pressure(data)
level0to1.assign_ucass_lut(data)
level0to1.bin_centre_dp_um(data)
level0to1.ucass_sample_volume(data)
level0to1.mass_concentration_umm3(data)

up_split = data.up_profile_mask
