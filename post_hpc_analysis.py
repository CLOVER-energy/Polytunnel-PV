#!/usr/bin/python3.10
########################################################################################
# post_hpc_analysis.py - Module for running analysis on HPC bypass-diode runs.         #
#                                                                                      #
# Author: Ben Winchester                                                               #
# Copyright: Ben Winchester, 2024                                                      #
# Date created: 12/12/2024                                                             #
# License: Open source                                                                 #
########################################################################################
"""
post_hpc_analysis.py - The post-HPC analysis module for Polytunnel-PV.

This module plots a bunch of things.

"""

import argparse
import json
import math
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import os
import seaborn as sns
import sys
import warnings

from datetime import datetime, timedelta
from collections import defaultdict
from matplotlib import pyplot as plt
from matplotlib import rc, rcParams
from typing import Any, Callable, Generator, Hashable, Match, Pattern

import numpy as np
import pandas as pd

from tqdm import tqdm
from tqdm.std import tqdm as tqdm_pbar

from src.polytunnelpv.__main__ import (
    _parse_cells,
    _parse_locations,
    _parse_polytunnel_curves,
    _parse_pv_modules,
    _parse_pv_system,
    _parse_scenarios,
    _parse_weather,
    _sanitise_time,
    _solar_angles_from_weather_row,
)
from src.polytunnelpv.__utils__ import NAME, VOLTAGE_RESOLUTION
from src.polytunnelpv.pv_module.bypass_diode import BypassDiode, BypassedCellString
from src.polytunnelpv.pv_module.pv_cell import (
    get_irradiance,
    PVCell,
    relabel_cell_electrical_parameters,
    ZERO_CELSIUS_OFFSET,
)
from src.polytunnelpv.pv_module.pv_module import (
    Curve,
    CurveType,
    CurvedPVModule,
    ModuleType,
    TYPE_TO_CURVE_MAPPING,
)
from src.polytunnelpv.pv_system import ModuleString, PVSystem
from src.polytunnelpv.scenario import Scenario

rc("font", **{"family": "sans-serif", "sans-serif": ["Arial"]})
# rc("figure", **{"figsize": (48 / 5, 32 / 5)})
rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42
sns.set_context("talk")
sns.set_style("ticks")

import warnings

warnings.filterwarnings("ignore")


def _parse_args(unparsed_args: list[Any]) -> argparse.Namespace:
    """
    Parse the CLI args.

    :param: unparsed_args
        The unparsed command-line arguments.

    :returns:
        The parsed CLI args.

    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--output-file",
        "-o",
        default=None,
        help="The name of the combined output file.",
        type=str,
    )

    return parser.parse_args(unparsed_args)


def main(args: list[Any]) -> None:
    """
    Main plotter function.

    :param: args
        The CLI args.

    """

    parsed_args = _parse_args(args)

    if (combined_output_filename := parsed_args.output_file) is None:
        raise argparse.ArgumentError(
            "--output-file", "The output file must be specified."
        )

    with open(combined_output_filename, "r", encoding="UTF-8") as output_file:
        combined_data = pd.read_csv(output_file, index_col=0)

    locations = _parse_locations()
    polytunnels = _parse_polytunnel_curves()
    user_defined_pv_cells = _parse_cells()
    pv_modules = _parse_pv_modules(
        {polytunnel.name: polytunnel for polytunnel in polytunnels},
        {entry[NAME]: entry for entry in user_defined_pv_cells},
    )
    pv_system = _parse_pv_system()
    scenarios = _parse_scenarios(
        {location.name: location for location in locations},
        {module.name: module for module in pv_modules},
    )

    # Remove ludicrous values
    combined_data.clip(lower=0, upper=5000, inplace=True)
    combined_data = combined_data[
        [entry for entry in combined_data if combined_data[entry].max(axis=0) < 5000]
    ]

    # Determine scenario information relevant for analysis.
    relevant_scenarios = {
        entry.name: entry for entry in scenarios if entry.name in combined_data.columns
    }
    bypass_diodes = {
        scenario.name: [
            entry
            for entry in scenario.pv_module.pv_cells_and_cell_strings
            if not isinstance(entry, PVCell)
        ]
        for scenario in relevant_scenarios.values()
    }
    num_bypass_diodes = {
        scenario_name: len(diodes) for scenario_name, diodes in bypass_diodes.items()
    }

    num_diodes_frame = pd.DataFrame(
        {"scenario": num_bypass_diodes.keys(), "num_diodes": num_bypass_diodes.values()}
    )
    num_diodes_frame.index = num_diodes_frame["scenario"]
    num_diodes_frame.pop("scenario")

    # Determine the impact of the number of bypass diodes, regardless of length
    diodes_to_frame: dict[int, pd.DataFrame] = {
        num_diodes: combined_data[
            [
                entry
                for entry in num_diodes_frame.index
                if num_diodes_frame.loc[entry, "num_diodes"] == num_diodes
            ]
        ]
        for num_diodes in set(num_bypass_diodes.values())
    }
    sns.set_palette(
        sns.color_palette(
            [
                "#8F1838",
                "#C11F33",
                "#EA1D2D",
                "#F36E24",
                "#F99D25",
                "#FDB714",
                "#00ACD7",
                "#007DBB",
                "#00558A",
                "#1A3668",
                "#48773C",
                "#40AE49",
            ]
        )
    )

    sns.set_palette(
        sns.cubehelix_palette(
            start=0.6, rot=-0.6, n_colors=len(set(diodes_to_frame.keys()))
        )
    )

    plt.figure(figsize=(48 / 5, 32 / 5))
    for num_diodes, sub_frame in reversed(diodes_to_frame.items()):
        sns.scatterplot(
            sub_frame.mean(axis=1),
            alpha=0.55,
            label=num_diodes,
            linewidth=0,
            marker="h",
            s=250,
        )

    norm = plt.Normalize(
        min(set(diodes_to_frame)) - 1,
        max(set(diodes_to_frame)) + 1,
    )
    scalar_mappable = plt.cm.ScalarMappable(
        cmap=mcolors.LinearSegmentedColormap.from_list(
            "Custom", sns.color_palette().as_hex(), len(set(diodes_to_frame)) + 1
        ),
        norm=norm,
    )

    colorbar = (axis := plt.gca()).figure.colorbar(
        scalar_mappable,
        ax=axis,
        label="Number of bypass diodes installed",
        pad=(_pad := 0.025),
    )

    plt.legend().remove()

    plt.xlabel("Hour of the day")
    plt.ylabel("System-wide power produced / kW")
    plt.savefig(
        f"hourly_power_scatter_{os.path.basename(combined_output_filename).replace('.csv', '')}.png",
        bbox_inches="tight",
        pad_inches=0.035,
        transparent=True,
    )
    plt.savefig(
        f"hourly_power_scatter_{os.path.basename(combined_output_filename).replace('.csv', '')}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )

    plt.show()

    diodes_to_power = {
        num_diodes: pd.DataFrame(
            {"power": frame.sum(axis=0).transpose().reset_index(drop=True)}
        )
        for num_diodes, frame in diodes_to_frame.items()
    }
    for num_diodes, frame in diodes_to_power.items():
        frame["num_diodes"] = num_diodes

    num_diodes_power_frame = pd.concat(diodes_to_power.values(), axis=0)

    # Bypass didoe number impact
    cat_plot = sns.catplot(
        num_diodes_power_frame,
        kind="boxen",
        x="num_diodes",
        y="power",
        height=9.75,
        aspect=(52 / 48),
        palette=sns.color_palette(),
        linewidth=0,
        linecolor=".7",
        line_kws=dict(linewidth=0, color="#cde"),
        showfliers=True,
    )

    plt.xlabel("Number of bypass diodes installed")
    plt.ylabel("System-wide power produced / kWh")
    # norm = plt.Normalize(
    #     min(set(diodes_to_frame)) - 1,
    #     max(set(diodes_to_frame)) + 1,
    # )
    # scalar_mappable = plt.cm.ScalarMappable(
    #     cmap=mcolors.LinearSegmentedColormap.from_list(
    #         "Custom", sns.color_palette().as_hex(), len(set(diodes_to_frame)) + 1
    #     ),
    #     norm=norm,
    # )

    # colorbar = (axis := plt.gca()).figure.colorbar(
    #     scalar_mappable,
    #     ax=axis,
    #     label="Number of bypass diodes installed",
    #     pad=(_pad := 0.025),
    # )

    plt.legend().remove()
    # plt.ylim(0, None)

    medians = num_diodes_power_frame.groupby(["num_diodes"])["power"].median()
    vertical_offset = (
        num_diodes_power_frame["power"].median() * 0.01
    )  # offset from median for display

    for xtick, xtick_position in zip(
        cat_plot.ax.get_xticklabels(), cat_plot.ax.get_xticks()
    ):
        if (xtick_value := int(xtick._text)) < plt.xlim()[0]:
            continue
        if xtick_value % 4 != 0:
            continue
        cat_plot.ax.text(
            xtick_position,
            medians[xtick_value] + vertical_offset,
            f"{medians[xtick_value]:.0f}",
            horizontalalignment="center",
            size="x-small",
            color=sns.color_palette().as_hex()[-1],
            # weight="semibold",
        )

    for ax in plt.gcf().axes:
        median = num_diodes_power_frame.groupby("num_diodes")["power"].median()
        for index, med in enumerate(median.values):
            ax.axhline(
                xmin=(x_pos := index / len(median.values)),
                y=med,
                xmax=x_pos + 1 / len(median.values),
                color="white",
            )

    plt.xticks(*[entry[::2] for entry in plt.xticks()])

    plt.savefig(
        f"number_of_diodes_{os.path.basename(combined_output_filename).replace('.csv', '')}.png",
        bbox_inches="tight",
        pad_inches=0.35,
        transparent=True,
    )
    plt.savefig(
        f"number_of_diodes_{os.path.basename(combined_output_filename).replace('.csv', '')}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )

    plt.show()

    cat_plot = sns.catplot(
        num_diodes_power_frame,
        kind="boxen",
        x="num_diodes",
        y="power",
        height=6.5,
        aspect=(48 / 32),
        palette=sns.color_palette(),
        linewidth=0,
        linecolor=".7",
        line_kws=dict(linewidth=0, color="#cde"),
        showfliers=True,
    )

    plt.xlabel("Number of bypass diodes installed")
    plt.ylabel("System-wide power produced / kWh")
    # norm = plt.Normalize(
    #     min(set(diodes_to_frame)) - 1,
    #     max(set(diodes_to_frame)) + 1,
    # )
    # scalar_mappable = plt.cm.ScalarMappable(
    #     cmap=mcolors.LinearSegmentedColormap.from_list(
    #         "Custom", sns.color_palette().as_hex(), len(set(diodes_to_frame)) + 1
    #     ),
    #     norm=norm,
    # )

    # colorbar = (axis := plt.gca()).figure.colorbar(
    #     scalar_mappable,
    #     ax=axis,
    #     label="Number of bypass diodes installed",
    #     pad=(_pad := 0.025),
    # )

    plt.legend().remove()

    plt.plot(
        (
            max_frame := num_diodes_power_frame.groupby(["num_diodes"])["power"]
            .max()
            .reset_index(drop=True)
        ).index,
        max_frame.values,
        color="grey",
        dashes=(2, 2),
    )

    plt.xticks(*[entry[::2] for entry in plt.xticks()])

    plt.savefig(
        f"number_of_diodes_with_max_{os.path.basename(combined_output_filename).replace('.csv', '')}.png",
        bbox_inches="tight",
        pad_inches=0.35,
        transparent=True,
    )
    plt.savefig(
        f"number_of_diodes_with_max_{os.path.basename(combined_output_filename).replace('.csv', '')}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )

    plt.show()

    # Plot staggered box plots
    plt.figure(figsize=(48 / 5, 32 / 5))
    melted_frames = []
    for num_diodes, sub_frame in reversed(diodes_to_frame.items()):
        melted_frame = sub_frame.reset_index().melt(
            id_vars="index", var_name="scenario"
        )
        melted_frame["bypass_diodes"] = num_diodes
        melted_frames.append(melted_frame)

    sns.boxplot(
        combined_melted_frame := pd.concat(melted_frames, axis=0),
        x="index",
        y="value",
        hue="bypass_diodes",
        dodge=True,
        palette=sns.color_palette(),
    )

    plt.xlim(5, 19)

    norm = plt.Normalize(
        min(set(diodes_to_frame)) - 1,
        max(set(diodes_to_frame)) + 1,
    )
    scalar_mappable = plt.cm.ScalarMappable(
        cmap=mcolors.LinearSegmentedColormap.from_list(
            "Custom", sns.color_palette().as_hex(), len(set(diodes_to_frame)) + 1
        ),
        norm=norm,
    )

    colorbar = (axis := plt.gca()).figure.colorbar(
        scalar_mappable,
        ax=axis,
        label="Number of bypass diodes installed",
        pad=(_pad := 0.025),
    )

    plt.legend().remove()
    plt.show()

    joint_plot = sns.jointplot(
        combined_melted_frame,
        x="index",
        y="value",
        alpha=0,
        hue="bypass_diodes",
        palette=sns.color_palette(),
        marginal_ticks=True,
        height=6.5,
        ratio=6,
    )

    sns.scatterplot(
        combined_melted_frame,
        x="index",
        y="value",
        hue="bypass_diodes",
        alpha=0.15,
        ax=joint_plot.ax_joint,
        # cbar=True,
        # cbar_kws={"label": "Irradiance / kWm$^{-2}$"},
        marker="D",
        palette=sns.color_palette(),
        s=150,
    )

    norm = plt.Normalize(
        min(set(diodes_to_frame)) - 0.5,
        max(set(diodes_to_frame)) + 0.5,
    )
    scalar_mappable = plt.cm.ScalarMappable(
        cmap=mcolors.LinearSegmentedColormap.from_list(
            "Custom", sns.color_palette().as_hex(), len(set(diodes_to_frame)) + 1
        ),
        norm=norm,
    )

    colorbar = (axis := joint_plot.ax_joint).figure.colorbar(
        scalar_mappable,
        ax=axis,
        label="Number of bypass diodes installed",
        pad=(_pad := 0.025),
    )

    plt.xlabel("Hour of the day")
    plt.ylabel("System-wide power produced / kW")

    joint_plot.ax_marg_x.remove()

    joint_plot.ax_joint.legend().remove()
    joint_plot.ax_joint.set_ylim(0, 265)
    joint_plot.ax_joint.set_xlim(4, 20)

    sns.despine(ax=joint_plot.ax_joint)
    sns.despine(ax=joint_plot.ax_marg_y)
    plt.subplots_adjust(
        # top=1.1,
        # bottom=0.1,
        # left=0.1,
        right=1.1
    )
    plt.savefig(
        f"bypass_diode_joint_plot_{os.path.basename(combined_output_filename).replace('.csv', '')}.png",
        bbox_inches="tight",
        pad_inches=0.35,
        transparent=True,
    )
    plt.savefig(
        f"bypass_diode_joint_plot_{os.path.basename(combined_output_filename).replace('.csv', '')}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    sorted_combined_data = combined_data[
        combined_data.sum(axis=0).sort_values(ascending=False).index
    ]
    sorted_combined_data.pop("hour")
    sorted_combined_data.columns = [
        f"Scenario {entry.split('n_')[1].split('_4344')[0].split('_4368')[0]}"
        for entry in sorted_combined_data.columns
    ]

    # Plot the top performing scenarios
    plt.figure(figsize=(48 / 5, 32 / 5))
    sns.lineplot(
        sliced_data := sorted_combined_data[
            reversed(sorted_combined_data.columns[: len(sns.color_palette())])
        ]
    )

    for index, column in enumerate(sliced_data.columns):
        plt.fill_between(
            sliced_data.index,
            sorted_combined_data["Scenario 21"],
            sliced_data[column],
            color=f"C{index}",
            alpha=0.1,
        )

    plt.xlabel("Hour of the day")
    plt.ylabel("System-wide power produced / kW")

    plt.xlim(4, 20)
    plt.ylim(0, None)

    plt.legend(ncols=2, title="Scenario ID", bbox_to_anchor=(1.1, 1.05))
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.savefig(
        f"top_scenarios_{os.path.basename(combined_output_filename).replace('.csv', '')}.png",
        bbox_inches="tight",
        pad_inches=0.35,
        transparent=True,
    )
    plt.savefig(
        f"top_scenarios_{os.path.basename(combined_output_filename).replace('.csv', '')}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )

    plt.show()

    plt.figure(figsize=(48 / 5, 32 / 5))

    for index, column in enumerate(sliced_data.columns):
        plt.fill_between(
            sliced_data.index,
            sorted_combined_data["Scenario 21"],
            sliced_data[column],
            color=f"C{index}",
            alpha=0.1,
        )

    plt.xlabel("Hour of the day")
    plt.ylabel("System-wide power produced / kW")

    plt.xlim(4, 20)
    plt.ylim(0, None)

    plt.legend(handles, labels, ncols=2, title="Scenario ID")

    plt.savefig(
        f"top_scenarios_fill_only_{os.path.basename(combined_output_filename).replace('.csv', '')}.png",
        bbox_inches="tight",
        pad_inches=0.35,
        transparent=True,
    )
    plt.savefig(
        f"top_scenarios_fill_only_{os.path.basename(combined_output_filename).replace('.csv', '')}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    # Plot the next set of scenarios
    plt.figure(figsize=(48 / 5, 32 / 5))
    sns.lineplot(
        sliced_data := sorted_combined_data[
            reversed(
                sorted_combined_data.columns[
                    len(sns.color_palette()) : 2 * len(sns.color_palette())
                ]
            )
        ]
    )

    for index, column in enumerate(sliced_data.columns):
        plt.fill_between(
            sliced_data.index,
            sorted_combined_data["Scenario 21"],
            sliced_data[column],
            color=f"C{index}",
            alpha=0.1,
        )

    plt.xlabel("Hour of the day")
    plt.ylabel("System-wide power produced / kW")

    plt.xlim(4, 20)
    plt.ylim(0, None)

    plt.legend(ncols=2, title="Scenario ID", bbox_to_anchor=(1.1, 1.05))
    handles, labels = plt.gca().get_legend_handles_labels()

    plt.savefig(
        f"second_scenarios_{os.path.basename(combined_output_filename).replace('.csv', '')}.png",
        bbox_inches="tight",
        pad_inches=0.35,
        transparent=True,
    )
    plt.savefig(
        f"second_scenarios_{os.path.basename(combined_output_filename).replace('.csv', '')}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    plt.figure(figsize=(48 / 5, 32 / 5))
    for index, column in enumerate(sliced_data.columns):
        plt.fill_between(
            sliced_data.index,
            sorted_combined_data["Scenario 21"],
            sliced_data[column],
            color=f"C{index}",
            alpha=0.1,
        )

    plt.xlabel("Hour of the day")
    plt.ylabel("System-wide power produced / kW")

    plt.xlim(4, 20)
    plt.ylim(0, None)

    plt.legend(handles, labels, ncols=2, title="Scenario ID")

    plt.savefig(
        f"second_scenarios_fill_only_{os.path.basename(combined_output_filename).replace('.csv', '')}.png",
        bbox_inches="tight",
        pad_inches=0.35,
        transparent=True,
    )
    plt.savefig(
        f"second_scenarios_fill_only_{os.path.basename(combined_output_filename).replace('.csv', '')}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    # Plot the next set of scenarios
    plt.figure(figsize=(48 / 5, 32 / 5))
    sns.lineplot(
        sliced_data := sorted_combined_data[
            reversed(
                sorted_combined_data.columns[
                    2 * len(sns.color_palette()) : 3 * len(sns.color_palette())
                ]
            )
        ]
    )

    for index, column in enumerate(sliced_data.columns):
        plt.fill_between(
            sliced_data.index,
            sorted_combined_data["Scenario 21"],
            sliced_data[column],
            color=f"C{index}",
            alpha=0.1,
        )

    plt.xlabel("Hour of the day")
    plt.ylabel("System-wide power produced / kW")

    plt.xlim(4, 20)
    plt.ylim(0, None)

    plt.legend(ncols=2, title="Scenario ID", bbox_to_anchor=(1.1, 1.05))
    handles, labels = plt.gca().get_legend_handles_labels()

    plt.savefig(
        f"third_scenarios_{os.path.basename(combined_output_filename).replace('.csv', '')}.png",
        bbox_inches="tight",
        pad_inches=0.35,
        transparent=True,
    )
    plt.savefig(
        f"third_scenarios_{os.path.basename(combined_output_filename).replace('.csv', '')}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    plt.figure(figsize=(48 / 5, 32 / 5))

    for index, column in enumerate(sliced_data.columns):
        plt.fill_between(
            sliced_data.index,
            sorted_combined_data["Scenario 21"],
            sliced_data[column],
            color=f"C{index}",
            alpha=0.1,
        )

    plt.xlabel("Hour of the day")
    plt.ylabel("System-wide power produced / kW")

    plt.xlim(4, 20)
    plt.ylim(0, None)

    plt.legend(handles, labels, ncols=2, title="Scenario ID")

    plt.savefig(
        f"third_scenarios_fill_only_{os.path.basename(combined_output_filename).replace('.csv', '')}.png",
        bbox_inches="tight",
        pad_inches=0.35,
        transparent=True,
    )
    plt.savefig(
        f"third_scenarios_fill_only_{os.path.basename(combined_output_filename).replace('.csv', '')}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.show()

    import pdb

    pdb.set_trace()

    for num_diodes, sub_frame in combined_frame.groupby(
        combined_frame.loc["num_diodes"]
    ):
        sns.catplot(sub_frame, kind="box")

    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
