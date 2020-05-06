from matplotlib import pyplot as plt
from matplotlib import font_manager as plt_fnt
from matplotlib import rcParams as mplParams
from matplotlib import lines
from matplotlib.legend import Legend
import matplotlib.ticker as plticker
from AirborneParticleAnalysis import common


plt.style.use("ggplot")
prop = plt_fnt.FontProperties(family=['serif'])
mplParams["font.family"] = prop.get_name()
mplParams['hatch.linewidth'] = 0.5
mplParams['mathtext.default'] = "regular"


def plot_rebin_1to1(data_1, data_2, regression_data, names):

    fig = plt.figure()
    fig.set_size_inches(common.cm_to_inch(8.3, 8.3))
    ax = fig.add_axes([0.2, 0.15, 0.75, 0.7])
    title_string = "Integrated dn/dlog(Dp)"
    ax.set_title(title_string, fontsize="small")
    ax.set_ylabel(names[1], fontsize="small")
    ax.set_xlabel(names[0], fontsize="small")

    x12 = regression_data[0]
    y12 = regression_data[1]
    r2 = regression_data[2]
    m = regression_data[5]

    marker_style = dict(linestyle='none', marker='x', markersize=5, fillstyle='none', color=(0, 0, 0), linewidth=0.7)

    ax.plot(data_1, data_2, **marker_style)

    marker_style["linestyle"] = '-.'
    marker_style["marker"] = ''
    ax.plot(x12, y12, **marker_style)
    h1 = lines.Line2D([], [], **marker_style)

    marker_style["linestyle"] = 'solid'
    ax.plot([0, x12[-1]], [0, x12[-1]], **marker_style)
    h2 = lines.Line2D([], [], **marker_style)

    ax.text(0.1, 0.9, "$r^{2} = $%f\n$m = $%f" % (r2, m),
            fontsize="small", ha='center', va='center', transform=ax.transAxes)

    ax.set_ylim(ymin=0)
    ax.set_xlim(xmin=0)

    leg = Legend(ax, [h1, h2], ["Regression Line", "y = x"], frameon=False, fontsize="small", loc=2)
    ax.add_artist(leg)

    plt.show()

    return
