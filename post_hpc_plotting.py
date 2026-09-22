#!/usr/bin/python3

import argparse
import datetime
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import os
import re
import sys

from matplotlib import rc
from typing import Any

import json
import numpy as np
import pandas as pd
import seaborn as sns

from tqdm import tqdm

# ALPHA:
#   The transparency to use for plotting.
ALPHA: float = 0.15

# DPI:
#   The DPI to use when saving figures.
DPI: int = 400

# INDEX:
#   Variable used for plotting to ensure non-overwriting behaviour when saving plots
INDEX: int = 4

# MM:
#   Conversion factor from mm to inches.
MM = 1 / 25.4

# MONTHS:
#   A list of names of months
MONTHS: list[str] = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]

PURPLE_GREEN_COLOURS: list[str] = [
    "#003D30",
    "#005745",
    "#00735C",
    "#009175",
    "#00AD8E",
    "#00CBA7",
    "#00EBC1",
    "#86FFDE",
    "#FFCCFE",
    "#FF92FD",
    "#FF3CFE",
    "#DA00FD",
    "#A700FC",
    "#8400CD",
    "#8400CD",
    "#65019F",
    "#450270",
]

# Plotting context
rc("font", **{"family": "sans-serif", "sans-serif": ["Arial"], "size": 7})
sns.set_context("paper", rc={"font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7})
sns.set_style("ticks")

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

from matplotlib import rcParams

rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42

plt.rcParams["font.size"] = 7


TOTEX_HEADER: str = "Other TOTEX"

rc("font", **{"family": "sans-serif", "sans-serif": ["Arial"]})
rcParams["text.usetex"] = False
sns.set_context("paper")
sns.set_style("ticks")

# Set custom color-blind colormap
colorblind_palette = sns.color_palette(
    [
        "#E04606",  # Orange
        "#F09F52",  # Pale orange
        "#52C0AD",  # Pale green
        "#006264",  # Green
        "#D8247C",  # Pink
        "#EDEDED",  # Pale pink
        "#E7DFBE",  # Pale yellow
        "#FBBB2C",  # Yellow
    ]
)

sns.set_palette(colorblind_palette)

thesis_palette = sns.color_palette(
    [
        "#E04606",
        "#F9A130",
        "#FCB919",  # SDG 7
        "#EFF2DD",
        "#27BFE6",  # SDG 6
        "#144E56",
        "#D8247C",  # Pink
        "#EDEDED",  # Pale pink
        "#E7DFBE",  # Pale yellow
    ]
)

categorical_colourblind_palette = sns.color_palette(
    [
        "#77AADD",
        "#99DDFF",
        "#44BB99",
        "#BBCC33",
        "#AAA000",
        "#EEDD88",
        "#EE8866",
        "#FFAABB",
        "#DDDDDD",
    ]
)

blend_palette = sns.blend_palette(["#36C7B8", "#423252"], n_colors=50)


sns.set_palette(
    blue_red_palette := sns.color_palette(
        [
            "#3455CC",
            "#557BF9",
            "#9FB0FC",
            "#D1C3A7",
            "#FBC412",
            "#FE8224",
            "#E03944",
        ]
    )
)

# Line of code for heatmap plotting
# cbar_axes = sns.heatmap(
#     pd.DataFrame(
#         {
#             int(hour): entry.mean()
#             for hour, entry in cigs_full_bypassing_frame.groupby("hour")
#         }
#     )
#     .transpose()
#     .set_index("hour"),
#     cmap="cividis",
#     cbar_kws={"label": "Average bypassing rate"},
# )
# (cbar_axis := cbar_axes.figure.axes[-1]).yaxis.label.set_size(7)
# cbar_axes.collections[0].colorbar.ax.tick_params(labelsize=7)
# plt.xlabel("Index of bypass-diode cell pair", fontsize=7)
# plt.ylabel("Hour", fontsize=7)
# (axis := plt.gca()).tick_params(axis="both", which="major", labelsize=7)
# plt.show()


def parse_args(args: list[Any]) -> argparse.Namespace:
    """
    Parse the CLI args.

    :param: args:
        The unparsed CLI args.

    :returns:
        The parsed args.

    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--cigs-output-file",
        "-cigs",
        type=str,
        default="12_aug_26_mpp_hourly_cigs.csv",
        help="The hourly output of the curved CIGS modules.",
    )
    parser.add_argument(
        "--m-si-output-file",
        "-msi",
        type=str,
        default="13_aug_m_si_pv_performance.csv",
        help="The hourly output of the flat m-Si modules.",
    )

    parser.add_argument(
        "--base-time",
        "-bt",
        type=str,
        default="2016-01-01 00:00",
        help="The base start time to use.",
    )

    parser.add_argument(
        "--array-size",
        type=int,
        default=56,
        help="The number of PV modules in the array.",
    )

    parser.add_argument(
        "--cigs-module-size",
        type=float,
        default=280,
        help="The capacity, in kWp, of a CIGS module.",
    )

    return parser.parse_args(args)


def main(args: list[Any]) -> None:
    """
    Main function.

    :param: args:
        The unparsed CLI args.

    """

    parsed_args = parse_args(args)

    with open(parsed_args.cigs_output_file, "r", encoding="UTF-8") as cigs_file:
        cigs_data = pd.read_csv(cigs_file, index_col=0)

    cigs_data.columns = pd.Index([_p_cigs_label := "P_cigs"])

    with open(parsed_args.m_si_output_file, "r", encoding="UTF-8") as msi_file:
        msi_data = pd.read_csv(msi_file, index_col=0)

    # Offset the CIGS data by the time.
    base_time = datetime.datetime.strptime(parsed_args.base_time, "%Y-%m-%d %H:%M")
    cigs_data["time"] = [
        base_time + datetime.timedelta(hours=entry) for entry in cigs_data.index
    ]

    # Set the m-Si data time correctly.
    try:
        msi_data["time"] = [
            datetime.datetime.strptime(entry, "%Y-%m-%d %H:%M")
            for entry in msi_data["time"]
        ]
    except ValueError:
        try:
            msi_data["time"] = [
                datetime.datetime.strptime(entry, "%Y-%m-%d %H:%M:%S")
                for entry in msi_data["time"]
            ]
        except ValueError as error:
            raise Exception("Could not process time information.") from error

    # Correct all units to kWp.
    cigs_data[_p_cigs_label] /= 1000

    ##########################
    # Array-wide power plots #
    ##########################

    sns.set_palette(categorical_colourblind_palette)
    plt.figure(figsize=(171 * MM, 120 * MM))
    sns.lineplot(cigs_data, x="time", y=_p_cigs_label, color="C0", label="PI-PV CIGS")
    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)
    (ax := plt.gca()).fill_between(
        cigs_data["time"],
        [0] * len(cigs_data),
        cigs_data[_p_cigs_label],
        alpha=ALPHA,
        color="C0",
    )
    ax.tick_params(axis="both", which="major", labelsize=7)

    plt.legend(fontsize=7)
    plt.savefig(
        f"cigs_output_power_10_yr_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    plt.figure(figsize=(83 * MM, 60 * MM))
    sns.lineplot(cigs_data, x="time", y=_p_cigs_label, color="C0", label="PI-PV CIGS")
    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)
    (ax := plt.gca()).fill_between(
        cigs_data["time"],
        [0] * len(cigs_data),
        cigs_data[_p_cigs_label],
        alpha=ALPHA,
        color="C0",
    )
    ax.tick_params(axis="both", which="major", labelsize=7)

    plt.legend(fontsize=7)
    plt.savefig(
        f"cigs_output_power_10_yr_small_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    sns.set_palette(categorical_colourblind_palette)
    plt.figure(figsize=(171 * MM, 120 * MM))
    sns.lineplot(cigs_data, x="time", y=_p_cigs_label, color="C0", label="PI-PV CIGS")
    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)
    (ax := plt.gca()).fill_between(
        cigs_data["time"],
        [0] * len(cigs_data),
        cigs_data[_p_cigs_label],
        alpha=ALPHA,
        color="C0",
    )
    ax.tick_params(axis="both", which="major", labelsize=7)
    plt.xlim(datetime.datetime(2024, 1, 1), datetime.datetime(2024, 12, 31, 23))

    plt.legend(fontsize=7)
    plt.savefig(
        f"cigs_output_power_2024_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    plt.figure(figsize=(83 * MM, 60 * MM))
    sns.lineplot(cigs_data, x="time", y=_p_cigs_label, color="C0", label="PI-PV CIGS")
    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)
    (ax := plt.gca()).fill_between(
        cigs_data["time"],
        [0] * len(cigs_data),
        cigs_data[_p_cigs_label],
        alpha=ALPHA,
        color="C0",
    )
    ax.tick_params(axis="both", which="major", labelsize=7)
    plt.xlim(datetime.datetime(2024, 1, 1), datetime.datetime(2024, 12, 31, 23))

    plt.legend(fontsize=7)
    plt.savefig(
        f"cigs_output_power_2024_small_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    sns.set_palette(categorical_colourblind_palette)
    plt.figure(figsize=(171 * MM, 120 * MM))
    sns.lineplot(cigs_data, x="time", y=_p_cigs_label, color="C0", label="PI-PV CIGS")
    (ax := plt.gca()).fill_between(
        cigs_data["time"],
        [0] * len(cigs_data),
        cigs_data[_p_cigs_label],
        alpha=ALPHA,
        color="C0",
    )
    ax.tick_params(axis="both", which="major", labelsize=7)
    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)
    plt.xlim(datetime.datetime(2024, 6, 1), datetime.datetime(2024, 6, 30, 23))

    plt.legend(fontsize=7)
    plt.savefig(
        f"cigs_output_power_june_2024_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    plt.figure(figsize=(83 * MM, 60 * MM))
    sns.lineplot(cigs_data, x="time", y=_p_cigs_label, color="C0", label="PI-PV CIGS")
    (ax := plt.gca()).fill_between(
        cigs_data["time"],
        [0] * len(cigs_data),
        cigs_data[_p_cigs_label],
        alpha=ALPHA,
        color="C0",
    )
    ax.tick_params(axis="both", which="major", labelsize=7)
    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)
    plt.xlim(datetime.datetime(2024, 6, 1), datetime.datetime(2024, 6, 30, 23))

    plt.legend(fontsize=7)
    plt.savefig(
        f"cigs_output_power_june_2024_small_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    sns.set_palette(categorical_colourblind_palette)
    plt.figure(figsize=(171 * MM, 120 * MM))
    sns.lineplot(cigs_data, x="time", y=_p_cigs_label, color="C0", label="PI-PV CIGS")
    (ax := plt.gca()).fill_between(
        cigs_data["time"],
        [0] * len(cigs_data),
        cigs_data[_p_cigs_label],
        alpha=ALPHA,
        color="C0",
    )
    ax.tick_params(axis="both", which="major", labelsize=7)
    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)
    plt.xlim(datetime.datetime(2025, 4, 1), datetime.datetime(2025, 4, 30, 23))

    plt.legend(fontsize=7)
    plt.savefig(
        f"cigs_output_power_apr_2025_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    plt.figure(figsize=(83 * MM, 60 * MM))
    sns.lineplot(cigs_data, x="time", y=_p_cigs_label, color="C0", label="PI-PV CIGS")
    (ax := plt.gca()).fill_between(
        cigs_data["time"],
        [0] * len(cigs_data),
        cigs_data[_p_cigs_label],
        alpha=ALPHA,
        color="C0",
    )
    ax.tick_params(axis="both", which="major", labelsize=7)
    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)
    plt.xlim(datetime.datetime(2025, 4, 1), datetime.datetime(2025, 4, 30, 23))

    plt.legend(fontsize=7)
    plt.savefig(
        f"cigs_output_power_apr_2025_small_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    ###############################
    # Array-wide power comparison #
    ###############################

    # Rescale the m-Si modules to match the output power of the array.
    capacity_column_headers = [
        re.match(r".*/P_max", entry).group(0)
        for entry in msi_data.columns
        if re.match(r".*/P_max", entry) is not None
    ]
    for column_name in capacity_column_headers:
        msi_data[column_name.replace("/P_max", "_array")] = (
            msi_data[column_name]
            * parsed_args.array_size
            * parsed_args.cigs_module_size
            / 1000
        )

    sns.set_palette(categorical_colourblind_palette)
    plt.figure(figsize=(171 * MM, 120 * MM))
    sns.lineplot(cigs_data, x="time", y=_p_cigs_label, color="C0", label="PI-PV CIGS")
    # (ax := plt.gca()).fill_between(
    #     cigs_data["time"],
    #     [0] * len(cigs_data),
    #     cigs_data[_p_cigs_label],
    #     alpha=ALPHA,
    #     color="C0",
    # )

    _index: int = 2
    for entry in msi_data.columns:
        if "array" not in entry:
            continue
        sns.lineplot(
            msi_data,
            x="time",
            y=entry,
            color=f"C{_index}",
            label=entry.split("P_")[1].split("_array")[0],
        )
        # ax.fill_between(msi_data["time"], [0] * len(msi_data), msi_data[entry], alpha=ALPHA, color=f"C{_index}")
        _index += 1

    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)

    plt.legend(fontsize=7)
    plt.gca().tick_params(axis="both", which="major", labelsize=7)

    plt.savefig(
        f"output_power_10_yr_all_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    plt.figure(figsize=(83 * MM, 60 * MM))
    sns.lineplot(cigs_data, x="time", y=_p_cigs_label, color="C0", label="PI-PV CIGS")
    # (ax := plt.gca()).fill_between(
    #     cigs_data["time"],
    #     [0] * len(cigs_data),
    #     cigs_data[_p_cigs_label],
    #     alpha=ALPHA,
    #     color="C0",
    # )

    _index: int = 2
    for entry in msi_data.columns:
        if "array" not in entry:
            continue
        sns.lineplot(
            msi_data,
            x="time",
            y=entry,
            color=f"C{_index}",
            label=entry.split("P_")[1].split("_array")[0],
        )
        # ax.fill_between(msi_data["time"], [0] * len(msi_data), msi_data[entry], alpha=ALPHA, color=f"C{_index}")
        _index += 1

    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)

    plt.gca().tick_params(axis="both", which="major", labelsize=7)
    plt.legend(fontsize=7)
    plt.savefig(
        f"output_power_10_yr_all_small_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    # Plot the output power of the individual types of m-Si modules.
    plotting_data = msi_data.set_index("time")[
        [entry for entry in msi_data.columns if "array" in entry]
    ]
    plt.figure(figsize=(171 * MM, 120 * MM))
    sns.lineplot(cigs_data, x="time", y=_p_cigs_label, color="C0", label="PI-PV CIGS")
    sns.lineplot(plotting_data.reset_index().set_index("time"))

    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)
    plt.gca().tick_params(axis="both", which="major", labelsize=7)

    plt.legend(fontsize=7)
    plt.savefig(
        f"output_power_10_yr_all_improved_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    plt.figure(figsize=(83 * MM, 60 * MM))
    sns.lineplot(cigs_data, x="time", y=_p_cigs_label, color="C0", label="PI-PV CIGS")
    sns.lineplot(plotting_data.reset_index().set_index("time"))

    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)
    plt.gca().tick_params(axis="both", which="major", labelsize=7)

    plt.legend(fontsize=7)
    plt.savefig(
        f"output_power_10_yr_all_improved_small_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    # Aggregate the output power of the modules.
    _fixed_tilt = pd.DataFrame(
        plotting_data[
            [entry for entry in plotting_data.columns if "fixed_tilt" in entry]
        ].mean(axis=1)
    )
    _fixed_tilt.columns = pd.Index(["Fixed-tilt"])
    _dual_axis = pd.DataFrame(
        plotting_data[
            [entry for entry in plotting_data.columns if "dual_axis" in entry]
        ].mean(axis=1)
    )
    _dual_axis.columns = pd.Index(["Dual-axis"])
    plotting_data = pd.merge(
        pd.merge(cigs_data.set_index("time"), _fixed_tilt, on="time"),
        _dual_axis,
        on="time",
    )
    plotting_data.columns = pd.Index(
        [entry.replace("P_cigs", "PI-PV CIGS") for entry in plotting_data.columns]
    )

    plt.figure(figsize=(171 * MM, 120 * MM))
    sns.lineplot(plotting_data)

    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)

    plt.legend(fontsize=7)
    plt.gca().tick_params(axis="both", which="major", labelsize=7)
    plt.savefig(
        f"output_power_10_yr_design_average_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    plt.figure(figsize=(83 * MM, 60 * MM))
    sns.lineplot(plotting_data)

    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)

    plt.legend(fontsize=7)
    plt.gca().tick_params(axis="both", which="major", labelsize=7)
    plt.savefig(
        f"output_power_10_yr_design_average_small_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    # Plot hourly
    plotting_data_hourly = plotting_data.copy()
    plotting_data_hourly["hour"] = [entry.hour for entry in plotting_data_hourly.index]
    plotting_data_hourly = pd.DataFrame(
        [entry.mean() for _, entry in plotting_data_hourly.groupby("hour")]
    ).set_index("hour")
    plotting_data_hourly.columns = pd.Index(
        [
            entry.replace("P_cigs", "PI-PV CIGS")
            for entry in plotting_data_hourly.columns
        ]
    )

    plt.figure(figsize=(171 * MM, 120 * MM))
    sns.lineplot(plotting_data_hourly)

    axis = plt.gca()
    for index, column in enumerate(plotting_data_hourly.columns):
        axis.fill_between(
            plotting_data_hourly.index,
            [0] * len(plotting_data_hourly),
            plotting_data_hourly[column],
            alpha=ALPHA,
            color=f"C{index}",
        )

    axis.tick_params(axis="both", which="major", labelsize=7)
    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)

    plt.legend(fontsize=7)
    plt.savefig(
        f"output_power_10_yr_hourly_average_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    plt.figure(figsize=(83 * MM, 60 * MM))
    sns.lineplot(plotting_data_hourly)

    axis = plt.gca()
    for index, column in enumerate(plotting_data_hourly.columns):
        axis.fill_between(
            plotting_data_hourly.index,
            [0] * len(plotting_data_hourly),
            plotting_data_hourly[column],
            alpha=ALPHA,
            color=f"C{index}",
        )

    axis.tick_params(axis="both", which="major", labelsize=7)
    plt.xlabel("Time", fontsize=7)
    plt.ylabel("Array power / kW", fontsize=7)

    plt.legend(fontsize=7)
    plt.savefig(
        f"output_power_10_yr_hourly_average_small_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    plotting_data_monthly = plotting_data.copy()
    plotting_data_monthly["month"] = [
        entry.month for entry in plotting_data_monthly.index
    ]
    plotting_data_monthly["hour"] = [
        entry.hour for entry in plotting_data_monthly.index
    ]
    plotting_data_monthly.columns = pd.Index(
        [
            entry.replace("P_cigs", "PI-PV CIGS")
            for entry in plotting_data_monthly.columns
        ]
    )
    monthly_means: dict[Any, pd.DataFrame] = {
        month: pd.DataFrame(
            [entry.mean() for _, entry in frame.groupby("hour")]
        ).set_index("hour")
        for month, frame in plotting_data_monthly.groupby("month")
    }
    for key, value in monthly_means.items():
        value.pop("month")
        monthly_means[key] = value

    sns.set_style("ticks")
    fig, axes = plt.subplots(4, 3, figsize=(171 * MM, 240 * MM))
    fig.subplots_adjust(hspace=0.35, wspace=0.25)
    flat_axes = (entry for entry in np.ndarray.flatten(axes))
    for month, frame in monthly_means.items():
        axis = next(flat_axes)
        sns.lineplot(frame, ax=axis)
        for index, column in enumerate(frame.columns):
            axis.fill_between(
                frame.index,
                [0] * len(frame),
                frame[column],
                alpha=ALPHA,
                color=f"C{index}",
            )
        axis.set_ylim(0, 12)
        axis.set_xlabel("")
        axis.set_ylabel("")
        axis.tick_params(axis="both", which="major", labelsize=7)
        axis.set_title(MONTHS[int(month) - 1], fontweight="bold", fontsize=7)
        if month < 12:
            axis.legend().remove()

    fig.text(0.5, 0.06, "Time", ha="center", fontsize=7)
    fig.text(
        0.04, 0.5, "Array power / kW", va="center", rotation="vertical", fontsize=7
    )

    handles, labels = axis.get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="upper left", bbox_to_anchor=(0.125, 0.875), fontsize=7
    )
    axis.legend().remove()

    plt.savefig(
        f"output_power_monthly_{INDEX}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.05,
    )
    plt.show()

    # Print out statistics for the compative output of the module designs over the
    # simulation period.
    print("Annual output kWh/kWp")
    print(
        (plotting_data.sum(axis=0) / (len(plotting_data) / 8760))
        / (parsed_args.array_size * (parsed_args.cigs_module_size / 1000))
    )
    print("Average daily output kWh/kWp")
    print(
        (plotting_data.sum(axis=0) * 24 / (len(plotting_data)))
        / (parsed_args.array_size * (parsed_args.cigs_module_size / 1000))
    )

    import pdb

    pdb.set_trace()


if __name__ == "__main__":
    main(sys.argv[1:])
