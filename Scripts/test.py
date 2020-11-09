from AirborneParticleAnalysis import level1to2
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
import numpy as np


if __name__ == "__main__":

    sam_data = level1to2.import_level1("C:\\Users\\jg17acv\\University of Hertfordshire\\AeroSAM - Documents\\Data"
                                       "\\2020\\20-09-28\\level_1\\FMITalon_PID-RM-003_20200928_10420556_00.pdat")
    t = sam_data.time[:3000]
    p = sam_data.press_hpa[:3000] * -1
    l, _ = find_peaks(np.squeeze(p), prominence=5, distance=10)
    plt.plot(p)
    plt.plot(l, p[l], "X")
    plt.show()
pass
