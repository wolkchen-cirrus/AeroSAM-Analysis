from AirborneParticleAnalysis import level1to2
import numpy as np


if __name__ == "__main__":
    data_file = "C:\\Users\\jg17acv\\University of Hertfordshire\\AeroSAM - Documents\\Data\\2020\\20-03-10\\level_1" \
                "\\CYI-FW-UCASS-X2_20200310_115555_00.pdat"
    out_file1 = "C:\\Users\\jg17acv\\University of Hertfordshire\\AeroSAM - Documents\\Data\\2020\\20-03-10\\level_2" \
                "\\115555_dvdlogdp_um3cm3.txt"
    out_file2 = "C:\\Users\\jg17acv\\University of Hertfordshire\\AeroSAM - Documents\\Data\\2020\\20-03-10\\level_2" \
                "\\115555_cint_dvdlogdp_um3cm3.txt"
    data = level1to2.import_level1(data_file)
    dn_slices = [200, 400, 600, 800, 1000]
    alt_list_up = []
    for dn_slice in dn_slices:
        alt_list_up.append(level1to2.fetch_row(altitude=dn_slice, level1_data=data, profile="Up"))
    alt_list_down = []
    for dn_slice in dn_slices:
        alt_list_down.append(level1to2.fetch_row(altitude=dn_slice, level1_data=data, profile="Down"))

    bins = [data.bin_centres_dp_um1, data.bin_centres_dp_um2]
    dn_dlogdp = [data.dv_dlogdp1.values(), data.dv_dlogdp2.values()]
    dn_dlogdp_arr = []
    for dn in dn_dlogdp:
        del dn[0:100]
        del dn[-100:-1]
        dn_arr = np.array([np.array(dni) for dni in dn])
        dn_dlogdp_arr.append(dn_arr)
    dn_p = []
    for dn in dn_dlogdp:
        dn_p.append(np.mean(dn, axis=0))
    with open(out_file2, "w+") as f:
        f.writelines("cintdv_ucass1 = %s\n" % (str(dn_p[0]).replace("\n", "")))
        f.writelines("cintdv_ucass2 = %s\n" % (str(dn_p[1]).replace("\n", "")))
        f.writelines("cint_bincent1 = %s\n" % (str(bins[0]).replace("\n", "")))
        f.writelines("cint_bincent2 = %s\n" % (str(bins[1]).replace("\n", "")))

    with open(out_file1, "w+") as f:
        for alt in alt_list_up:
            rows_for_mean = level1to2.fetch_row_tolerance(altitude=alt, level1_data=data, profile="Up")
            up_mean1, up_sigma1 = level1to2.mean_dn_dlogdp(data, rows_for_mean, ucass_number=1)
            up_mean2, up_sigma2 = level1to2.mean_dn_dlogdp(data, rows_for_mean, ucass_number=2)
            f.writelines("up_mean_ucass1_%s = %s\n" % (str(alt), str(up_mean1).replace("\n", "")))
            f.writelines("up_sigma_ucass1_%s = %s\n" % (str(alt), str(up_sigma1).replace("\n", "")))
            f.writelines("up_mean_ucass2_%s = %s\n" % (str(alt), str(up_mean2).replace("\n", "")))
            f.writelines("up_sigma_ucass2_%s = %s\n" % (str(alt), str(up_sigma2).replace("\n", "")))
        for alt in alt_list_down:
            rows_for_mean = level1to2.fetch_row_tolerance(altitude=alt, level1_data=data, profile="Down")
            down_mean1, down_sigma1 = level1to2.mean_dn_dlogdp(data, rows_for_mean, ucass_number=1)
            down_mean2, down_sigma2 = level1to2.mean_dn_dlogdp(data, rows_for_mean, ucass_number=2)
            f.writelines("down_mean_ucass1_%s = %s\n" % (str(alt), str(down_mean1).replace("\n", "")))
            f.writelines("down_sigma_ucass1_%s = %s\n" % (str(alt), str(down_sigma1).replace("\n", "")))
            f.writelines("down_mean_ucass2_%s = %s\n" % (str(alt), str(down_mean2).replace("\n", "")))
            f.writelines("down_sigma_ucass2_%s = %s\n" % (str(alt), str(down_sigma2).replace("\n", "")))
