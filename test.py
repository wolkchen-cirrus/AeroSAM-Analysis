import importer
import level0to1


data = importer.SUAData()
level0to1.split_by_pressure(data)

up_split = data.up_profile_mask
