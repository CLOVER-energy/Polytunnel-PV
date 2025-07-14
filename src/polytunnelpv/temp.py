import seaborn as sns
import matplotlib.pyplot as plt

import math

import matplotlib.colors as mcolors

sns.set_palette(
    "cividis",
    n_colors=(
        num_deviations := math.ceil(max(heat_pump_data["abs_deviation"]))
        - math.floor(min(heat_pump_data["abs_deviation"]))
    ),
)

plt.figure(figsize=(48 / 5, 32 / 5))
sns.scatterplot(
    heat_pump_data,
    x="heat_output",
    y="cop",
    hue="abs_deviation",
    palette="cividis",
    marker="D",
    s=250,
    alpha=0.5,
)
axis = plt.gca()

norm = plt.Normalize(
    min(heat_pump_data["abs_deviation"]),
    max(heat_pump_data["abs_deviation"]),
)
scalar_mappable = plt.cm.ScalarMappable(
    cmap=mcolors.LinearSegmentedColormap.from_list(
        "Custom", sns.color_palette().as_hex(), num_deviations
    ),
    norm=norm,
)

colorbar = (axis := plt.gca()).figure.colorbar(
    scalar_mappable,
    ax=axis,
    label="Absolute percentage deviation from model / %",
    pad=(_pad := 0.025),
)
plt.legend().remove()

plt.xlabel("Heat output at A7/W55 / kW")
plt.ylabel("COP")

plt.savefig(
    "heat_pump_cop_deviation_plot.pdf", format="pdf", bbox_inches="tight", pad_inches=0
)

plt.show()


plt.figure(figsize=(48 / 5, 32 / 5))

sns.set_palette(
    cmap := sns.cubehelix_palette(
        start=-0.4,
        rot=-0.6,
        n_colors=(
            num_deviations := math.ceil(10 * max(heat_pump_data["cop"]))
            - math.floor(10 * min(heat_pump_data["cop"]))
        ),
    )
)

sns.scatterplot(
    heat_pump_data,
    x="heat_output",
    hue="cop",
    y="abs_deviation",
    palette=cmap,
    marker="D",
    s=250,
    alpha=0.5,
)
axis = plt.gca()

norm = plt.Normalize(
    min(heat_pump_data["cop"]),
    max(heat_pump_data["cop"]),
)
scalar_mappable = plt.cm.ScalarMappable(
    cmap=mcolors.LinearSegmentedColormap.from_list(
        "Custom", sns.color_palette().as_hex(), num_deviations
    ),
    norm=norm,
)

colorbar = (axis := plt.gca()).figure.colorbar(
    scalar_mappable,
    ax=axis,
    label="COP",
    pad=(_pad := 0.025),
)
plt.legend().remove()

plt.xlabel("Heat output at A7/W55 / kW")
plt.ylabel("Absolute percentage deviation from model / %")

plt.axhline(4.78, dashes=(2, 1), color="grey", label="Mean absolute deviation")
plt.axhline(6, dashes=(2, 1, 4, 1), color="orange", label="Reported model undertainty")
plt.fill_between((xlim := [0, 27]), [0.0, 0.0], [0.2, 0.2], color="grey", alpha=0.2)
plt.fill_between(
    xlim,
    [9.36, 9.36],
    [20, 20],
    color="grey",
    alpha=0.2,
    label="Outside 1 standard deviation",
)
plt.xlim(*xlim)
plt.ylim(0, 18)

handles, labels = axis.get_legend_handles_labels()
plt.legend(handles[-3:], labels[-3:])

plt.savefig(
    "heat_pump_deviation_cop_plot.pdf", format="pdf", bbox_inches="tight", pad_inches=0
)

plt.show()
