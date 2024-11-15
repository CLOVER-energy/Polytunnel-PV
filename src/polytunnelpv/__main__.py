#!/usr/bin/python3.10
########################################################################################
# __main__.py - Main module for Polytunnel-PV.                                         #
#                                                                                      #
# Author: Ben Winchester                                                               #
# Copyright: Ben Winchester, 2024                                                      #
# Date created: 21/02/2024                                                             #
# License: Open source                                                                 #
########################################################################################
"""
__main__.py - The main module for Polytunnel-PV.

Polytunnel-PV simulates the performance of curved photovoltaic modules for polytunnel
and greenhouse applications. This main module provides a command-line interface
entrypoint for executing the model.

"""

__version__ = "1.0.0a2"


# Use AGG when running on HPC - uncomment below if running on HPC
# matplotlib.use('Agg')
import argparse
import enum
import functools
import json
import math
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import multiprocessing
import os
import pvlib
import re
import seaborn as sns
import sys
import time
import warnings
import yaml

from datetime import datetime, timedelta
from collections import defaultdict
from collections.abc import Iterator
from joblib import Parallel, delayed
from matplotlib import pyplot as plt
from matplotlib import rc, rcParams
from multiprocessing.dummy import Pool as ThreadPool
from typing import Any, Callable, Generator, Hashable, Match, Pattern

import numpy as np
import pandas as pd

from tqdm import tqdm
from tqdm.std import tqdm as tqdm_pbar

from .__utils__ import NAME, VOLTAGE_RESOLUTION
from .pv_module.bypass_diode import BypassDiode, BypassedCellString
from .pv_module.pv_cell import (
    get_irradiance,
    PVCell,
    relabel_cell_electrical_parameters,
    ZERO_CELSIUS_OFFSET,
)
from .pv_module.pv_module import (
    Curve,
    CurveType,
    CurvedPVModule,
    ModuleType,
    TYPE_TO_CURVE_MAPPING,
)
from .pv_system import ModuleString, PVSystem
from .scenario import Scenario

rc("font", **{"family": "sans-serif", "sans-serif": ["Arial"]})
# rc("figure", **{"figsize": (48 / 5, 32 / 5)})
rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42
sns.set_context("notebook")
sns.set_style("ticks")

# BYPASS_DIODES:
#   Keyword for the bypass-diode parameters.
BYPASS_DIODES: str = "bypass_diodes"

# CELL_ELECTRICAL_PARAMETERS:
#   Keyword for the electrical parameters of a PV cell.
CELL_ELECTRICAL_PARAMETERS: str = "cell_electrical_parameters"

# CELL_TYPE:
#   Keyword for the name of the cell-type to use.
CELL_TYPE: str = "cell_type"

# CURRENT_SAMPLING_RATE:
#    Rate used to sample current data for ease of calculation.
CURRENT_SAMPLING_RATE: int = 20

# DONE:
#   The message to display when a task was successful.
DONE: str = "[   DONE   ]"

# FAILED:
#   The message to display when a task was successful.
FAILED: str = "[  FAILED  ]"

# FILE_ENCODING:
#   The encoding to use when opening and closing files.
FILE_ENCODING: str = "UTF-8"

# HOUR:
#   Column header to use when saving hour information.
HOUR: str = "Hour"

# INPUT_DATA_DIRECTORY:
#   The name of the input-data directory.
INPUT_DATA_DIRECTORY: str = "input_data"

# INTERNAL_AXIS_PADDING_FACTOR:
#   Factor used to provide space between limits and axis extremes.
INTERNAL_AXIS_PADDING_FACTOR: float = 0.05

# IRRADIANCE_DIFFUSE:
#   Keyword for diffuse irradiance.
IRRADIANCE_DIFFUSE: str = "irradiance_diffuse"

# IRRADIANCE_DIRECT:
#   Keyword for direct irradiance.
IRRADIANCE_DIRECT: str = "irradiance_direct"

# IRRADIANCE_DIRECT_NORMAL:
#   Keyword for the DNI.
IRRADIANCE_DIRECT_NORMAL: str = "irradiance_direct_normal"

# IRRADIANCE_GLOBAL_HORIZONTAL:
#   Keyword for the GHI.
IRRADIANCE_GLOBAL_HORIZONTAL: str = "irradiance_global_horizontal"

# LOCAL_TIME:
#   Column header for local-time column.
LOCAL_TIME: str = "local_time"

# LOCATIONS_FILENAME:
#   The filename for the locations file.
LOCATIONS_FILENAME: str = "locations.yaml"

# MODULE_STRINGS:
#   Keyword for the module strings.
MODULE_STRINGS: str = "module_strings"

# N_MODULES:
#   The number of modules in the string.
N_MODULES: str = "n_modules"

# OUTPUT_DIRECTORY:
#   The name of the outputs directory to use.
OUTPUT_DIRECTORY: str = "outputs"

# POLYTUNNEL_CURVE@
#   Keyword used for parsing the information about the curve on which a solar panel
# sits.
POLYTUNNEL_CURVE: str = "polytunnel_curve"

# POLYTUNNEL_HEADER_STRING:
#   Header string for the polytunnel-PV code.
POLYTUNNEL_HEADER_STRING: str = """
                                   #         #        #
                                   ###       ##     ###
                            ##      ##################      #
                             ####  #########################
                                ##########################
                          #######################################
                              ##############################
                              ##############################
                                               #############
                                                    ########
                              ################           #
                          #######           #####
 ###   ###              #####     ########     #####
 ##### ###             ###      #####  #######    ####
 ### #####            ###          ########         #####
 ###   ###           ###                   #################
                     ##               #######             ####
        ##            ###          #####         #######    #####
        #####          ###       ####        ######      ###   ####
         #######        ####    ###                #########      ####
          #####           #######                ###                 ####
          ## ####          #####                       #####################
              ####           ###                   ######              #######
                ####           ###             ######                        #####
                  ###           ####         #####                             #####
                   ####           ###       ###                                   ###
                     ###           ####    ###                                     ####
                      ####           ###  ##                                        ###
                        ###           ######                                         ###
                         ####           ###                                          ###
                           ####
                            ####
                              ###
                               ####
                                 ####
                                  ###

{version_line}
                                     Polytunnel-PV
    An open-source modelling framework for the simulation and optimisation of curved
                      photovoltaic panels in agricultural contexts

                             For more information, contact
                  Benedict Winchester (benedict.winchester@gmail.com)
"""

# POLYTUNNELS_FILENAME:
#   The name of the polytunnels file.
POLYTUNNELS_FILENAME: str = "polytunnels.yaml"

# PV_CELL:
#   Keyword for parsing cell-id information.
PV_CELL: str = "cell_id"

# PV_CELLS_FILENAME:
#   The name of the PV-cells file.
PV_CELLS_FILENAME: str = "pv_cells.yaml"

# PV_SYSTEM_FILENAME:
#   The name of the PV-system file.
PV_SYSTEM_FILENAME: str = "pv_system.yaml"

# PVLIB_DATABASE_NAME:
#   The database to use when fetching data on PV cells from pvilb.
PVLIB_DATABASE_NAME: str = "CECmod"

# PV_MODULES_FILENAME:
#   The name of the PV-modules file.
PV_MODULES_FILENAME: str = "pv_modules.yaml"

# SCENARIOS_FILENAME:
#   The name of the scenarios file.
SCENARIOS_FILENAME: str = "scenarios.yaml"

# SKIPPED:
#   The message to display when a task was successful.
SKIPPED: str = "[  SKIPPED ]"

# SOLAR_AZIMUTH:
#   Keyword for solar azimuth.
SOLAR_AZIMUTH: str = "azimuth"

# SOLAR_ZENITH:
#   Keyword for apparent zenith.
SOLAR_ZENITH: str = "apparent_zenith"

# TEMPERATURE:
#   Column header for temperature column.
TEMPERATURE: str = "temperature"

# TIME:
#   Column header for time.
TIME: str = "time"

# TYPE:
#   Keyword used to determine the module type of the PV.
TYPE: str = "type"

# Version regex:
#   Regex used to extract the main version number.
VERSION_REGEX: Pattern[str] = re.compile(r"(?P<number>\d\.\d\.\d)([\.](?P<post>.*))?")

# WEATHER_DATA_DIRECTORY:
#   The directory where weather data should be found.
WEATHER_DATA_DIRECTORY: str = "weather_data"

# WEATHER_DATA_REGEX:
#   Regex used for parsing location names from weather data.
WEATHER_DATA_REGEX: Pattern[str] = re.compile(r"ninja_pv_(?P<location_name>.*)\.csv")

# WIND_DATA_REGEX:
#   Regex used for parsing location names from weather data.
WIND_DATA_REGEX: Pattern[str] = re.compile(r"ninja_wind_(?P<location_name>.*)\.csv")

# WIND_SPEED:
#   Column header for wind-speed column.
WIND_SPEED: str = "wind_speed"

# WEATHER_FILE_WITH_SOLAR:
#   The name of the weather file with the solar data.
WEATHER_FILE_WITH_SOLAR: str = os.path.join(
    "auto_generated", "weather_with_solar_{location_name}.csv"
)


class ArgumentError(Exception):
    """Raised when incorrect command-line arguments are passed in."""

    def __init__(self, msg: str) -> None:
        """Instantiate with the message `msg`."""

        super().__init__(f"Incorrect CLI arguments: {msg}")


class Dashes(Iterator):
    """
    Contains the dash information for plotting.

    .. attribute:: dashes
        The dashes to include.

    """

    def __init__(self) -> None:
        """
        Instantiate the dashes.

        """

        self.dash_length: int = 0
        self.space: int = 1
        super().__init__()

    def __next__(self) -> Generator[tuple[int, int], None, None]:
        """
        Iterate through and yield the dash information.

        :yields: The dash information as a `tuple`.

        """

        # Move on the dash length
        self.dash_length += 1

        # Move on the space length if needed
        if self.dash_length - self.space > 3:
            self.dash_length = 1
            self.space += 1

        return (self.dash_length, self.space)


class OperatingMode(enum.Enum):
    """
    Used to determine what operation to carry out.

    - HOURLY_MPP:
        Carry out a computation for the hourly MPP across the module.

    - HOURLY_MPP_HPC:
        Carry out a computation for the hourly MPP across the module when operating on a
        high-performance computer (HPC), *i.e.*, in parallel operation.

    - IRRADIANCE_HEATMAPS:
        Generate irradiance and temperature heatmaps.

    - IRRADIANCE_LINEPLOTS:
        Generate irradiance and temperature lineplots.

    - IV_CURVES:
        Generate IV curves for the module.

    - VALIDATION:
        Validate the model against inputted data.

    """

    HOURLY_MPP: str = "hourly_mpp"
    HOURLY_MPP_HPC: str = "hourly_mpp_hpc"
    HOURLY_MPP_MACHINE_LEARNING_TRAING: str = "hourly_mpp_machine_learning"
    IRRADIANCE_HEATMAPS: str = "irradiance_heatmaps"
    IRRADIANCE_LINEPLOTS: str = "irradiance_lineplots"
    IV_CURVES: str = "iv_curves"
    VALIDATION: str = "validation"


def _parse_args(unparsed_args: list[str]) -> argparse.Namespace:
    """
    Parse the command-line arguments.

    :param: **unparsed_args:**
        The unparsed command-line arguments.

    """

    parser = argparse.ArgumentParser()

    # Iteration length:
    #   The length of the iteration to run, in hours.
    parser.add_argument(
        "--iteration-length",
        "-i",
        default=24,
        type=int,
        help="The length of the iteration to run.",
    )

    # Scenario:
    #   The name of the scenario to use.
    parser.add_argument(
        "--scenario", "-s", type=str, help="The name of the scenario to use."
    )

    # Start-day index:
    #   The index of the start day to use.
    parser.add_argument(
        "--start-day-index",
        "-st",
        default=0,
        type=int,
        help="The start-day index to use.",
    )

    # Timestamps file:
    #   The name of the timestamps file to use for determining when to simulate.
    parser.add_argument(
        "--timestamps-file",
        "-t",
        type=str,
        help="The name of the timestamps file to use for simualtions.",
        default=None,
    )

    # Operating mode:
    #   The operation to carry out.
    parser.add_argument(
        "--operating-mode",
        choices=[entry.value for entry in OperatingMode],
        help="Specify the operating mode.",
    )

    # Weather-data arguments.
    weather_arguments = parser.add_argument_group("weather-data arguments")
    # Regenerate:
    #   Used to re-calculate the solar weather information.
    weather_arguments.add_argument(
        "--regenerate",
        "-rw",
        action="store_true",
        default=False,
        help="Recalculate the weather and solar information.",
    )

    return parser.parse_args(unparsed_args)


def _parse_cells() -> list[dict[str, float]]:
    """
    Parses the PV cell parameters based on the input file.

    PV cells can either be fetched from a database of `pvlib` information or specified
    by the user. If the user is specifying the cell parameters, then this is done in the
    input file.

    """

    with open(
        os.path.join(INPUT_DATA_DIRECTORY, PV_CELLS_FILENAME),
        "r",
        encoding=FILE_ENCODING,
    ) as f:
        return yaml.safe_load(f)


def _parse_locations() -> list[pvlib.location.Location]:
    """Parses the locations based on the input file."""

    with open(
        os.path.join(INPUT_DATA_DIRECTORY, LOCATIONS_FILENAME),
        "r",
        encoding=FILE_ENCODING,
    ) as f:
        locations_data = yaml.safe_load(f)

    try:
        return [pvlib.location.Location(**entry) for entry in locations_data]
    except KeyError:
        raise KeyError("Not all location information present in locations file.")


def _parse_polytunnel_curves() -> list[Curve]:
    """Parse the polytunnel curves from the files."""

    with open(
        os.path.join(INPUT_DATA_DIRECTORY, POLYTUNNELS_FILENAME),
        "r",
        encoding=FILE_ENCODING,
    ) as f:
        polytunnels_data = yaml.safe_load(f)

    try:
        return [
            TYPE_TO_CURVE_MAPPING[CurveType(polytunnel_entry.pop(TYPE))](  # type: ignore [misc]
                **polytunnel_entry
            )
            for polytunnel_entry in polytunnels_data
        ]
    except KeyError:
        raise KeyError(
            f"Missing type entry with key '{TYPE}' for polytunnel curve."
        ) from None


def _parse_pv_modules(
    polytunnels: dict[str, Curve], user_defined_pv_cells: dict[str, dict[str, float]]
) -> list[CurvedPVModule]:
    """
    Parse the curved PV module information from the files.

    :param: **polytunnels:**
        A mapping between polytunnel names and instances.

    :param: **user_defined_pv_cells:**
            A mapping between PV cell name and pv cell information.

    :returns:
        The parsed PV modules as a list.

    """

    with open(
        os.path.join(INPUT_DATA_DIRECTORY, PV_MODULES_FILENAME),
        "r",
        encoding=FILE_ENCODING,
    ) as f:
        pv_module_data = yaml.safe_load(f)

    def _construct_pv_module(pv_module_entry) -> CurvedPVModule:
        """
        Constructs a PV module.

        :param: pv_module_entry
            The entry in the module file.

        """

        try:
            constructor = CurvedPVModule.constructor_from_module_type(
                ModuleType(pv_module_entry.pop(TYPE))
            )
        except KeyError:
            raise KeyError(
                f"Missing type entry with key '{TYPE}' for PV module."
            ) from None

        pv_module_entry[POLYTUNNEL_CURVE] = polytunnels[
            pv_module_entry[POLYTUNNEL_CURVE]
        ]

        # Attempt to determine the PV-cell parameters
        try:
            cell_type_name: str = pv_module_entry.pop(CELL_TYPE)
        except KeyError:
            raise KeyError(
                "Need to specify name of PV-cell material in pv-module file."
            ) from None

        # Try first to find the cell definition in a user-defined scope.
        try:
            cell_electrical_parameters = user_defined_pv_cells[cell_type_name]
        except KeyError:
            # Look within the pvlib database for the cell name.
            try:
                cell_electrical_parameters = (
                    pvlib.pvsystem.retrieve_sam(PVLIB_DATABASE_NAME)
                    .loc[:, cell_type_name]
                    .to_dict()
                )
            except KeyError:
                raise Exception(
                    f"Could not find cell name {cell_type_name} within either local or "
                    "pvlib-imported scope."
                ) from None

        # Create bypass diodes based on the information provided.
        try:
            pv_module_entry[BYPASS_DIODES] = (
                [BypassDiode(**entry) for entry in pv_module_entry[BYPASS_DIODES]]
                if BYPASS_DIODES in pv_module_entry
                else []
            )
        except KeyError as caught_error:
            raise Exception(
                "Missing bypass diode information for pv module."
            ) from caught_error

        # Map the parameters to those electrical parameters that the cell is expecting.
        pv_module_entry[CELL_ELECTRICAL_PARAMETERS] = (
            relabel_cell_electrical_parameters(cell_electrical_parameters)
        )

        return constructor(**pv_module_entry)

    return [_construct_pv_module(pv_module_entry) for pv_module_entry in pv_module_data]


def _parse_pv_system() -> PVSystem:
    """
    Parse the PV-system information from the input files.

    Outputs:
        The PV-system information.

    """

    with open(
        os.path.join(INPUT_DATA_DIRECTORY, PV_SYSTEM_FILENAME),
        "r",
        encoding=FILE_ENCODING,
    ) as f:
        pv_system_data = yaml.safe_load(f)

    try:
        return PVSystem(
            [ModuleString(entry[N_MODULES]) for entry in pv_system_data[MODULE_STRINGS]]
        )
    except KeyError:
        raise KeyError("Missing information in pv-system input file.") from None


def _parse_scenarios(
    locations: dict[str, pvlib.location.Location], pv_modules: dict[str, CurvedPVModule]
) -> list[Scenario]:
    """
    Parse the scenario information.

    :param: **locations:**
        The `list` of locations to use.

    :praam: **pv_modules:**
            The `list` of PVModules that can be installed at each location.

    :returns:
        - **scenarios:**
            The `list` of scenarios to run.

    """

    with open(
        os.path.join(INPUT_DATA_DIRECTORY, SCENARIOS_FILENAME),
        "r",
        encoding=FILE_ENCODING,
    ) as f:
        scenarios_data = yaml.safe_load(f)

    return [
        Scenario.from_scenarios_file(entry, locations, pv_modules)
        for entry in scenarios_data
    ]


def _parse_weather() -> dict[str, pd.DataFrame]:
    """
    Parse the downloaded solar data that's in the weather data directory.

    :returns: A `dict` mapping the location name, as a `str`, to a
        :class:`pandas.DataFrame` instance containing the solar and wind data.

    """

    location_name_to_data_map: dict[str, pd.DataFrame] = {}

    # Parse the solar data
    for filename in os.listdir(WEATHER_DATA_DIRECTORY):
        # Skip the file if it's not in the expected format.
        try:
            location_name = WEATHER_DATA_REGEX.match(filename).group("location_name")  # type: ignore [union-attr]
        except AttributeError:
            continue

        with open(
            os.path.join(WEATHER_DATA_DIRECTORY, filename), "r", encoding=FILE_ENCODING
        ) as f:
            location_name_to_data_map[location_name] = pd.read_csv(f, comment="#")

    for filename in os.listdir(WEATHER_DATA_DIRECTORY):
        # Skip the file if it's not in the expected format.
        try:
            location_name = WIND_DATA_REGEX.match(filename).group("location_name")  # type: ignore [union-attr]
        except AttributeError:
            continue

        with open(
            os.path.join(WEATHER_DATA_DIRECTORY, filename), "r", encoding=FILE_ENCODING
        ) as f:
            wind_data = pd.read_csv(f, comment="#")

        location_name_to_data_map[location_name][WIND_SPEED] = wind_data[WIND_SPEED]

    return location_name_to_data_map


def _sanitise_time(hours: int, formatter: str, initial_time: datetime) -> str:
    """
    Sanitise the time.

    :param: hours
        The hours to include.

    :param: formatter
        The string formatter to use.

    """

    return (initial_time + timedelta(hours=hours)).strftime(formatter)


def _solar_angles_from_weather_row(
    row: tuple[Hashable, pd.Series], location: pvlib.location.Location
) -> pd.DataFrame:
    """
    Use a row from the weather data to comopute the solar angles.

    :param: **row:**
        The row in the weather-data frame.

    :returns:
        The solar angle-frame at this time.

    """

    return location.get_solarposition(  # type: ignore [no-any-return]
        row[1][TIME], temperature=row[1][TEMPERATURE]
    )


def process_single_mpp_calculation(
    time_of_day: int,
    *,
    irradiance_frame: pd.DataFrame,
    locations_to_weather_and_solar_map: dict[
        pvlib.location.Location, list[dict[str, Any]]
    ],
    pbar: tqdm_pbar,
    pv_system: PVSystem,
    scenario: Scenario,
) -> tuple[int, float, str]:
    try:
        max_irradiance = np.max(
            irradiance_frame.set_index("hour").iloc[time_of_day][1:]
        )
        if max_irradiance == 0:
            return None, None, None

        hour = time_of_day % 24
        date = datetime(2023, 1, 1) + timedelta(hours=time_of_day)
        date_str = date.strftime("%d_%b_%Y")

        # Create a mapping between cell and power output
        cell_to_power_map: dict[BypassedCellString | PVCell, np.ndarray] = {}
        cell_to_voltage_map: dict[BypassedCellString | PVCell, np.ndarray] = {}
        cell_to_bypass_map: dict[BypassedCellString | PVCell, list[bool]] = {}

        current_series = np.linspace(
            0,
            1.1
            * max(
                [
                    pv_cell.short_circuit_current
                    for pv_cell in scenario.pv_module.pv_cells
                ]
            ),
            VOLTAGE_RESOLUTION,
        )

        individual_power_extreme = 0
        for pv_cell in tqdm(
            scenario.pv_module.pv_cells_and_cell_strings,
            desc="IV calculation",
            leave=False,
        ):
            (
                current_series,
                power_series,
                voltage_series,
                pv_cells_bypassed,
            ) = pv_cell.calculate_iv_curve(
                (
                    weather_map := locations_to_weather_and_solar_map[
                        scenario.location
                    ][time_of_day]
                )[TEMPERATURE],
                1000
                * irradiance_frame.set_index("hour")
                .iloc[time_of_day][1:]
                .reset_index(drop=True),
                weather_map[WIND_SPEED],
                current_series=current_series,
            )
            cell_to_power_map[pv_cell] = power_series
            cell_to_voltage_map[pv_cell] = voltage_series
            cell_to_bypass_map[pv_cell] = pv_cells_bypassed
            individual_power_extreme = max(
                individual_power_extreme, max(abs(power_series))
            )

        combined_power_series = sum(cell_to_power_map.values())
        combined_voltage_series = sum(cell_to_voltage_map.values())
        combined_power_series = pv_system.combine_powers(combined_power_series)
        combined_voltage_series = pv_system.combine_voltages(combined_voltage_series)
        combined_current_series = pv_system.combine_currents(current_series)

        maximum_power_index = pd.Series(combined_power_series).idxmax()
        mpp_current = combined_current_series[maximum_power_index]
        mpp_power = combined_power_series[maximum_power_index]

        pbar.update(1)

        return hour, mpp_power, date_str
    except Exception as e:
        print(f"Error processing time_of_day {time_of_day}: {str(e)}")
        raise
        # return None, None, None


def process_single_mpp_calculation_without_pbar(
    time_of_day: int,
    *,
    irradiance_frame: pd.DataFrame,
    locations_to_weather_and_solar_map: dict[
        pvlib.location.Location, list[dict[str, Any]]
    ],
    pv_system: PVSystem,
    scenario: Scenario,
) -> tuple[
    int,
    float,
    dict[BypassedCellString | PVCell, bool],
    dict[BypassedCellString | PVCell, float],
]:
    """
    Process the MPP calculation for a single hour of the simulation.

    :param: time_of_day
        The time of day to simulate.

    :param: irradiance_frame
        The irradiance data, as a :class:`pd.DataFrame`.

    :param: locations_to_weather_and_solar_map
        A `dict` mapping the location, as a :class:`pvlib.location.Location`, to the
        weather and solar data.

    :param: pv_system
        The :class:`PVSystem` instance representing the solar polytunnel system.

    :param: scenario
        The current :class:`Scenario` being modelled.

    :returns: The hour of the day,

    :returns: The MPP power achieved,

    :returns: Which bypassed cell strings were bypassing,

    :returns: and the MPP power achieved by each cell/string.

    """

    def _print_success() -> None:
        """Print a message when a computation was successful."""

        print(f"Hour {time_of_day} processed successfully.")

    try:
        hour = time_of_day % 24
        max_irradiance = np.max(
            irradiance_frame.set_index("hour").iloc[time_of_day][1:]
        )
        if max_irradiance == 0:
            _print_success()
            return (
                hour,
                0,
                {
                    cell_or_string: 0
                    for cell_or_string in scenario.pv_module.pv_cells_and_cell_strings
                },
                {
                    cell_or_string: 0
                    for cell_or_string in scenario.pv_module.pv_cells_and_cell_strings
                },
            )

        # Create a mapping between cell and power output
        cell_to_current_map: dict[BypassedCellString | PVCell, np.ndarray] = {}
        cell_to_power_map: dict[BypassedCellString | PVCell, np.ndarray] = {}
        cell_to_voltage_map: dict[BypassedCellString | PVCell, np.ndarray] = {}
        cell_to_bypass_map: dict[BypassedCellString | PVCell, list[bool]] = {}

        current_series = np.linspace(
            0,
            1.1
            * max(
                [
                    pv_cell.short_circuit_current
                    for pv_cell in scenario.pv_module.pv_cells
                ]
            ),
            VOLTAGE_RESOLUTION,
        )

        individual_power_extreme = 0
        for pv_cell in tqdm(
            scenario.pv_module.pv_cells_and_cell_strings,
            desc="IV calculation",
            leave=False,
        ):
            (
                current_series,
                power_series,
                voltage_series,
                pv_cells_bypassed,
            ) = pv_cell.calculate_iv_curve(
                (
                    weather_map := locations_to_weather_and_solar_map[
                        scenario.location
                    ][time_of_day]
                )[TEMPERATURE],
                1000
                * irradiance_frame.set_index("hour")
                .iloc[time_of_day]
                .reset_index(drop=True),
                weather_map[WIND_SPEED],
                current_series=current_series,
            )
            cell_to_current_map[pv_cell] = current_series
            cell_to_power_map[pv_cell] = power_series
            cell_to_voltage_map[pv_cell] = voltage_series
            cell_to_bypass_map[pv_cell] = pv_cells_bypassed
            individual_power_extreme = max(
                individual_power_extreme, max(abs(power_series))
            )

        # Compute the total power produced
        cell_voltage_interpreters = {
            pv_cell: lambda i, pv_cell=pv_cell: np.interp(
                i,
                list(reversed(cell_to_current_map[pv_cell])),
                list(reversed(cell_to_voltage_map[pv_cell])),
                period=np.max(cell_to_current_map[pv_cell])
                - np.min(cell_to_current_map[pv_cell]),
            )
            for pv_cell in scenario.pv_module.pv_cells_and_cell_strings
        }

        cellwise_voltage = {
            pv_cell: [
                interpreter(current)
                for current in tqdm(
                    current_series[::CURRENT_SAMPLING_RATE],
                    desc="Interpolation calculation",
                    leave=False,
                )
            ]
            for pv_cell, interpreter in tqdm(
                cell_voltage_interpreters.items(),
                desc="Cell-wise calculation",
                leave=False,
            )
        }
        module_voltage = [sum(sublist) for sublist in zip(*cellwise_voltage.values())]

        module_power = module_voltage * current_series[::CURRENT_SAMPLING_RATE]
        mpp_index: int = list(module_power).index(np.max(module_power))

        # Determine the power through the module at MPP
        mpp_power: float = module_power[mpp_index]

        # Determine the MPP power produced by each cell
        cellwise_mpp: dict[PVCell, float] = {
            pv_cell: voltage_interpreter(
                (mpp_current := current_series[::CURRENT_SAMPLING_RATE][mpp_index])
            )
            * mpp_current
            for pv_cell, voltage_interpreter in cell_voltage_interpreters.items()
        }

        def _find_nearest_index(input_list: list, target_value: float):
            """
            Determine the index of the nearest element to the target entry.

            :param: input_list
                The list which the entry should be looked up in.

            :param: target_value
                The value which to lookup.

            :returns: The index of the closest element.

            """

            return min(
                range(len(input_list)),
                key=lambda index: abs(input_list[index] - target_value),
            )

        # Determine which cell strings were bypassing
        bypassing_cell_strings: dict[PVCell, bool] = {
            pv_cell: bypass_map[
                _find_nearest_index(
                    cell_to_voltage_map[pv_cell],
                    cell_voltage_interpreters[pv_cell](mpp_current),
                )
            ]
            for pv_cell, bypass_map in cell_to_bypass_map.items()
        }

        _print_success()

        return hour, mpp_power, bypassing_cell_strings, cellwise_mpp

    except Exception as e:
        print(f"Error processing time_of_day {time_of_day}: {str(e)}")
        raise
        # return None, None, None


def plot_irradiance_with_marginal_means(
    data_to_plot: pd.DataFrame,
    start_index: int = 0,
    *,
    figname: str = "output_figure",
    fig_format: str = "pdf",
    heatmap_vmax: float | None = None,
    initial_date: datetime = datetime(2015, 1, 1),
    irradiance_bar_vmax: float | None = None,
    irradiance_scatter_vmax: float | None = None,
    irradiance_scatter_vmin: float | None = None,
    show_figure: bool = True,
) -> None:
    """
    Plot the irradiance heatmap with marginal means.

    :param: **data_to_plot:**
        The data to plot.

    :param: **start_index:**
        The index to start plotting from.

    :param: **figname:**
        The name to use when saving the figure.

    :param: **fig_format:**
        The file extension format to use when saving the high-resolution output figures.

    :param: **heatmap_vmax:**
        The maximum plotting value to use for the irradiance heatmap.

    :param: **irradiance_bar_vmax:**
        The maximum plotting value to use for the irradiance bars.

    :param: **irradiance_scatter_vmax:**
        The maximum plotting value to use for the irradiance scatter.

    :param: **irradiance_scatter_vmin:**
        The minimum plotting value to use for the irradiance scatter.

    """

    sns.set_style("ticks")
    frame_slice = (
        data_to_plot.fillna(0).iloc[start_index : start_index + 24].set_index("hour")
    )

    with tqdm(desc="Plotting", leave=False, total=3, unit="steps") as pbar:
        # Create a dummy plot and surrounding axes.
        joint_plot_grid = sns.jointplot(
            data=frame_slice,
            kind="hist",
            bins=(len(frame_slice.columns), 24),
            marginal_ticks=True,
            ratio=4,
        )
        joint_plot_grid.ax_marg_y.cla()
        joint_plot_grid.ax_marg_x.cla()
        pbar.update(1)

        # Generate the central heatmap.
        sns.heatmap(
            frame_slice,
            ax=joint_plot_grid.ax_joint,
            cmap=sns.blend_palette(
                ["#144E56", "teal", "#94B49F", "orange", "#E04606"], as_cmap=True
            ),
            vmin=0,
            vmax=heatmap_vmax,
            cbar=True,
            cbar_kws={"label": "Irradiance / kWm$^{-2}$"},
        )
        pbar.update(1)

        # Plot the irradiance bars on the RHS of the plot.
        joint_plot_grid.ax_marg_y.barh(
            np.arange(0.5, 24.5), frame_slice.mean(axis=1).to_numpy(), color="#E04606"
        )
        joint_plot_grid.ax_marg_y.set_xlim(0, math.ceil(10 * irradiance_bar_vmax) / 10)

        # Plot the scatter at the top of the plot.
        joint_plot_grid.ax_marg_x.scatter(
            np.arange(0.5, len(frame_slice.columns) + 0.5),
            (cellwise_means := frame_slice.mean(axis=0).to_numpy()),
            color="#144E56",
            marker="D",
            s=60,
            alpha=0.7,
        )
        joint_plot_grid.ax_marg_x.set_ylim(
            irradiance_scatter_vmin, irradiance_scatter_vmax
        )

        joint_plot_grid.ax_joint.set_xlabel("Mean angle from polytunnel axis")
        joint_plot_grid.ax_joint.set_ylabel("Hour of the day")
        joint_plot_grid.ax_joint.set_title(
            (initial_date + timedelta(hours=start_index)).strftime("%d/%m/%Y"),
            fontweight="bold",
        )
        joint_plot_grid.ax_joint.legend().remove()

        # Remove ticks from axes
        joint_plot_grid.ax_marg_x.tick_params(axis="x", bottom=False, labelbottom=False)
        joint_plot_grid.ax_marg_x.set_xlabel("Average irradiance / kWm$^{-2}$")
        joint_plot_grid.ax_marg_y.tick_params(axis="y", left=False, labelleft=False)
        # joint_plot_grid.ax_marg_y.set_xlabel("Average irradiance / kWm$^{-2}$")
        # remove ticks showing the heights of the histograms
        # joint_plot_grid.ax_marg_x.tick_params(axis='y', left=False, labelleft=False)
        # joint_plot_grid.ax_marg_y.tick_params(axis='frame_slice', bottom=False, labelbottom=False)

        joint_plot_grid.ax_marg_x.grid("on")
        # joint_plot_grid.ax_marg_x.set_ylabel("Average irradiance")  #

        joint_plot_grid.figure.tight_layout()

        # Adjust the plot such that the colourbar appears to the RHS.
        # Solution from JohanC (https://stackoverflow.com/users/12046409/johanc)
        # Obtained with permission from stackoverflow:
        # https://stackoverflow.com/questions/60845764/how-to-add-a-colorbar-to-the-side-of-a-kde-jointplot
        #
        plt.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)
        pos_joint_ax = joint_plot_grid.ax_joint.get_position()
        pos_marg_x_ax = joint_plot_grid.ax_marg_x.get_position()
        joint_plot_grid.ax_joint.set_position(
            [pos_joint_ax.x0, pos_joint_ax.y0, pos_marg_x_ax.width, pos_joint_ax.height]
        )
        joint_plot_grid.figure.axes[-1].set_position(
            [0.83, pos_joint_ax.y0, 0.07, pos_joint_ax.height]
        )

        # Save the figure in high- and low-resolution formats.
        plt.savefig(
            f"{figname}.{fig_format}",
            format=fig_format,
            bbox_inches="tight",
            pad_inches=0,
        )

        if show_figure:
            plt.show()

        pbar.update(1)


def plot_temperature_with_marginal_means(
    data_to_plot: pd.DataFrame,
    start_index: int = 0,
    *,
    figname: str = "output_figure",
    fig_format: str = "pdf",
    heatmap_vmax: float | None = None,
    initial_date: datetime = datetime(2015, 1, 1),
    locations_to_weather_and_solar_map: dict[
        pvlib.location.Location, list[dict[str, Any]]
    ],
    scenario: Scenario,
    show_figure: bool = True,
    temperature_bar_vmax: float | None = None,
    temperature_scatter_vmax: float | None = None,
    temperature_scatter_vmin: float | None = 0,
) -> None:
    """
    Plot the temperature heatmap with marginal means.

    :param: **data_to_plot:**
        The data to plot.

    :param: **start_index:**
        The index to start plotting from.

    :param: **figname:**
        The name to use when saving the figure.

    :param: **fig_format:**
        The file extension format to use when saving the high-resolution output figures.

    :param: **heatmap_vmax:**
        The maximum plotting value to use for the irradiance heatmap.

    :param: **initial_date:**
        The initial date for the run.

    :param: **locations_to_weather_and_solar_map:**
        A map between the location name and the weather and solar data.

    :param: **scenario:**
        The :class:`Scenario` governing the run.

    :param: **temperature_bar_vmax:**
        The maximum plotting value to use for the irradiance bars.

    :param: **temperature_scatter_vmax:**
        The maximum plotting value to use for the irradiance scatter.

    :param: **temperature_scatter_vmin:**
        The minimum plotting value to use for the irradiance scatter.

    """

    sns.set_style("ticks")
    frame_slice = (
        data_to_plot.fillna(0).iloc[start_index : start_index + 24].set_index("hour")
    )

    cell_to_temperature_profile_map: dict[float | int, list[float]] = {
        pv_cell.cell_id: [
            pv_cell.average_cell_temperature(
                locations_to_weather_and_solar_map[scenario.location][time_of_day][
                    TEMPERATURE
                ]
                + ZERO_CELSIUS_OFFSET,
                1000 * frame_slice.iloc[time_of_day].iloc[pv_cell.cell_id],
                0,
            )
            - ZERO_CELSIUS_OFFSET
            for time_of_day in tqdm(
                range(24), desc="Daily computation", leave=False, unit="hour"
            )
        ]
        for pv_cell in tqdm(
            scenario.pv_module.pv_cells,
            desc="Cell-wise computation",
            leave=False,
            unit="cell",
        )
    }
    temperature_frame = pd.DataFrame(cell_to_temperature_profile_map)
    temperature_frame.columns = frame_slice.columns

    if temperature_scatter_vmax is None:
        temperature_scatter_vmax = 1.25 * (
            max(
                temperature_frame.mean(axis=0).max(),
                0,
            )
        )

    with tqdm(desc="Plotting", leave=False, total=3, unit="steps") as pbar:
        # Create a dummy plot and surrounding axes.
        joint_plot_grid = sns.jointplot(
            data=temperature_frame,
            kind="hist",
            bins=(len(temperature_frame.columns), 24),
            marginal_ticks=True,
        )
        joint_plot_grid.ax_marg_y.cla()
        joint_plot_grid.ax_marg_x.cla()
        pbar.update(1)

        # Generate the central heatmap.
        sns.heatmap(
            temperature_frame,
            ax=joint_plot_grid.ax_joint,
            cmap="coolwarm",
            vmin=0,
            vmax=heatmap_vmax,
            cbar=True,
            cbar_kws={"label": "Temperature / $^\circ$C"},
        )
        pbar.update(1)

        # Plot the irradiance bars on the RHS of the plot.
        joint_plot_grid.ax_marg_y.barh(
            np.arange(0.5, 24.5), temperature_frame.mean(axis=1).to_numpy(), color="red"
        )
        joint_plot_grid.ax_marg_y.set_xlim(0, temperature_bar_vmax)

        # Plot the scatter at the top of the plot.
        joint_plot_grid.ax_marg_x.scatter(
            np.arange(0.5, len(temperature_frame.columns) + 0.5),
            (cellwise_means := temperature_frame.mean(axis=0).to_numpy()),
            color="red",
            marker="D",
            s=60,
            alpha=0.7,
        )
        joint_plot_grid.ax_marg_x.set_ylim(
            temperature_scatter_vmin, temperature_scatter_vmax
        )

        joint_plot_grid.ax_joint.set_xlabel("Mean angle from polytunnel axis")
        joint_plot_grid.ax_joint.set_ylabel("Hour of the day")
        joint_plot_grid.ax_joint.set_title(
            (initial_date + timedelta(hours=start_index)).strftime("%d/%m/%Y"),
            fontweight="bold",
        )
        joint_plot_grid.ax_joint.legend().remove()

        # Remove ticks from axes
        joint_plot_grid.ax_marg_x.tick_params(axis="x", bottom=False, labelbottom=False)
        joint_plot_grid.ax_marg_x.set_xlabel("Average temperature / $^\circ$C")
        joint_plot_grid.ax_marg_y.tick_params(axis="y", left=False, labelleft=False)
        # joint_plot_grid.ax_marg_y.set_xlabel("Average irradiance / kWm$^{-2}$")
        # remove ticks showing the heights of the histograms
        # joint_plot_grid.ax_marg_x.tick_params(axis='y', left=False, labelleft=False)
        # joint_plot_grid.ax_marg_y.tick_params(axis='frame_slice', bottom=False, labelbottom=False)

        joint_plot_grid.ax_marg_x.grid("on")
        # joint_plot_grid.ax_marg_x.set_ylabel("Average irradiance")  #

        joint_plot_grid.figure.tight_layout()

        # Adjust the plot such that the colourbar appears to the RHS.
        # Solution from JohanC (https://stackoverflow.com/users/12046409/johanc)
        # Obtained with permission from stackoverflow:
        # https://stackoverflow.com/questions/60845764/how-to-add-a-colorbar-to-the-side-of-a-kde-jointplot
        #
        plt.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)
        pos_joint_ax = joint_plot_grid.ax_joint.get_position()
        pos_marg_x_ax = joint_plot_grid.ax_marg_x.get_position()
        joint_plot_grid.ax_joint.set_position(
            [pos_joint_ax.x0, pos_joint_ax.y0, pos_marg_x_ax.width, pos_joint_ax.height]
        )
        joint_plot_grid.figure.axes[-1].set_position(
            [0.83, pos_joint_ax.y0, 0.07, pos_joint_ax.height]
        )

        # Save the figure in high- and low-resolution formats.
        plt.savefig(
            f"{figname}.{fig_format}",
            format=fig_format,
            bbox_inches="tight",
            pad_inches=0,
        )

        if show_figure:
            plt.show()

        pbar.update(1)


def save_daily_results(date_str, data, scenario_name):
    with pd.ExcelWriter(output_file, engine="openpyxl", mode="a") as writer:
        # df.to_excel(writer, sheet_name=date_str, index=False)
        # Ensure at least one sheet is present and visible
        workbook = writer.book
        placeholder_sheet = workbook.create_sheet(title="Sheet1")

        for date_str, data in daily_data.items():
            df = pd.DataFrame(data, columns=["Hour", "MPP (W)"])
            total_mpp = df["MPP (W)"].sum()
            df.loc["Total"] = ["Total", total_mpp]  # Adding the total MPP at the end
            df.to_excel(writer, sheet_name=date_str, index=False)

        # Remove the placeholder sheet if it was not used
        if "Sheet1" in workbook.sheetnames and len(workbook.sheetnames) > 1:
            del workbook["Sheet1"]


def main(unparsed_arguments) -> None:
    """
    Main method for Polytunnel-PV.

    :param: unparsed_arguments
        Un-parsed command-line arguments.

    """

    # Snippet taken with permission from CLOVER-energy/CLOVER
    # >>>
    version_match: Match[str] | None = VERSION_REGEX.match(__version__)
    version_number: str = (
        version_match.group("number") if version_match is not None else __version__
    )
    version_string = f"Version {version_number}"
    print(
        POLYTUNNEL_HEADER_STRING.format(
            version_line=(
                " " * (44 - math.ceil(len(version_string) / 2))
                + version_string
                + " " * (44 - math.floor(len(version_string) / 2))
            )
        )
    )
    # <<< end of reproduced snippted

    # Matplotlib setup
    rc("font", **{"family": "sans-serif", "sans-serif": ["Arial"]})
    sns.set_context("notebook")
    sns.set_style("ticks")

    # Parse the command-line arguments.
    parsed_args = _parse_args(unparsed_arguments)

    # Parse all of the input files
    print(
        (this_string := "Parsing input files")
        + "." * (88 - (len(this_string) + len(DONE)))
        + " ",
        end="",
    )
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

    # Parse the weather data.
    # NOTE: When integrated this as a Python package, this line should be suppressable
    # by weather data being passed in.
    weather_data = _parse_weather()
    print(DONE)

    # Map locations to weather data.
    print(
        (this_string := "Computing weather data")
        + "." * (88 - (len(this_string) + len(DONE)))
        + " ",
        end="",
    )
    locations_with_weather: dict[pvlib.location.Location, pd.DataFrame] = {
        location: weather_data[location.name]
        for location in locations
        if location.name in weather_data
    }

    # Map of locations to weather data with solar angles.
    locations_with_weather_and_solar: dict[pvlib.location.Location, pd.DataFrame] = {}

    for location, weather_frame in locations_with_weather.items():
        # Open the existing combined file with solar data if it exists
        if (
            os.path.isfile(
                (
                    weather_with_solar_filename := WEATHER_FILE_WITH_SOLAR.format(
                        location_name=location.name
                    )
                )
            )
            and not parsed_args.regenerate
        ):
            with open(weather_with_solar_filename, "r", encoding=FILE_ENCODING) as f:
                locations_with_weather_and_solar[location] = pd.read_csv(f, index_col=0)
        else:
            # Compute the solar-position information using the get-irradiance function.
            # with Pool(8) as worker_pool:
            #     solar_position_map = worker_pool.map(_solar_angles_from_weather_row, weather_frame.to_dict().items())
            # this_map = map(functools.partial(_solar_angles_from_weather_row, location=location), weather_frame.iterrows())
            solar_frame = pd.concat(
                [
                    _solar_angles_from_weather_row(row, location)
                    for row in tqdm(
                        weather_frame.iterrows(),
                        desc=location.name.capitalize(),
                        total=len(weather_frame),
                    )
                ]
            )

            # Compute the GHI from this information.
            weather_frame[IRRADIANCE_GLOBAL_HORIZONTAL] = (
                weather_frame[IRRADIANCE_DIRECT] + weather_frame[IRRADIANCE_DIFFUSE]
            )

            # Compute the DNI from this information.
            weather_frame[IRRADIANCE_DIRECT_NORMAL] = pvlib.irradiance.dni(
                weather_frame[IRRADIANCE_GLOBAL_HORIZONTAL].reset_index(drop=True),
                pd.Series(weather_frame[IRRADIANCE_DIFFUSE]).reset_index(drop=True),
                pd.Series(solar_frame[SOLAR_ZENITH]).reset_index(drop=True),
            )

            # Generate and save the combined frame.
            locations_with_weather_and_solar[location] = pd.concat(
                [solar_frame.reset_index(drop=True), weather_frame], axis=1
            )
            with open(weather_with_solar_filename, "w", encoding=FILE_ENCODING) as f:
                locations_with_weather_and_solar[location].to_csv(f)  # type: ignore[arg-type]

    # Convert each entry within the values of the mapping to a `dict` for speed.
    locations_to_weather_and_solar_map: dict[
        pvlib.location.Location, list[dict[str, Any]]
    ] = {
        location: entry.to_dict("records")
        for location, entry in locations_with_weather_and_solar.items()
    }
    print(DONE)

    # Compute the irradiance on each panel for each location.
    print(
        (this_string := "Calculting cell-wise irradiances")
        + "." * (88 - (len(this_string) + len(DONE)))
        + " ",
        end="",
    )
    # Load cell-wise irradiance data if already computed for the scenario being modelled
    # Otherwise, compute the irradiance data
    try:
        modelling_scenario = {scenario.name: scenario for scenario in scenarios}[
            parsed_args.scenario
        ]
    except KeyError:
        if parsed_args.scenario is None:
            raise ArgumentError(
                "Scenario must be specified on the command-line."
            ) from None
        raise KeyError(
            f"Scenario {parsed_args.scenario.name} not found in scenarios file. Valid scenarios: {', '.join([s.name for s in scenarios])}"
        ) from None

    cellwise_irradiance_frames: list[tuple[Scenario, pd.DataFrame]] = []

    if (
        not os.path.isfile(
            (
                cellwise_filename := os.path.join(
                    "auto_generated",
                    f"{modelling_scenario.name}_cellwise_irradiance.csv",
                )
            )
        )
        or parsed_args.regenerate
    ):
        irradiances_list = [
            {
                entry[TIME]: get_irradiance(
                    pv_cell,
                    entry[IRRADIANCE_DIFFUSE],
                    entry[IRRADIANCE_GLOBAL_HORIZONTAL],
                    entry[SOLAR_AZIMUTH],
                    entry[SOLAR_ZENITH],
                    direct_normal_irradiance=entry[IRRADIANCE_DIRECT_NORMAL],
                )
                for entry in tqdm(
                    locations_to_weather_and_solar_map[modelling_scenario.location][
                        :8760
                    ],
                    desc="Performing hourly calculation",
                    leave=False,
                )
            }
            for pv_cell in tqdm(
                modelling_scenario.pv_module.pv_cells,
                desc="Cell-wise irradiance",
                leave=False,
            )
        ]

        hourly_frames: list[pd.DataFrame] = []
        for cell_id, hourly_irradiances in enumerate(irradiances_list):
            hourly_frame = pd.DataFrame.from_dict(hourly_irradiances, orient="index")
            hourly_frame.columns = pd.Index([cell_id])
            hourly_frames.append(hourly_frame)

        combined_frame = pd.concat(hourly_frames, axis=1)
        combined_frame.sort_index(axis=1, inplace=True, ascending=False)

        combined_frame["hour"] = [
            int(entry.split(" ")[1].split(":")[0])
            for entry in hourly_frame.index.to_series()
        ]

        # Use the cell angle from the axis as the column header
        combined_frame.columns = [
            f"{'-' if pv_cell.cell_id / len(modelling_scenario.pv_module.pv_cells) < 0.5 else ''}{str(round(pv_cell.tilt, 3))}"
            for pv_cell in modelling_scenario.pv_module.pv_cells
        ] + ["hour"]

        with open(
            cellwise_filename,
            "w",
        ) as cellwise_file:
            combined_frame.to_csv(cellwise_file, index=None)

        cellwise_irradiance_frames.append((modelling_scenario, combined_frame))

        print(DONE)

    else:
        for scenario in tqdm(
            [scenario for scenario in scenarios if scenario == modelling_scenario],
            desc="Loading irradiance data",
            leave=True,
        ):
            with open(
                os.path.join(
                    "auto_generated", f"{scenario.name}_cellwise_irradiance.csv"
                ),
                "r",
            ) as f:
                combined_frame = pd.read_csv(f, index_col=None)

            # combined_frame.columns = pd.Index(
            #     [(date_and_time := "date_and_time")] + list(combined_frame.columns)[1:]
            # )
            cellwise_irradiance_frames.append((scenario, combined_frame.copy()))

        print(DONE)

    # Defragment the frame
    cellwise_irradiance_frames = [
        (entry[0], entry[1].copy()) for entry in cellwise_irradiance_frames
    ]

    try:
        irradiance_frame = [
            entry[1]
            for entry in cellwise_irradiance_frames
            if entry[0] == modelling_scenario
        ][0]
    except IndexError:
        raise Exception("Internal error occurred.") from None

    sns.set_palette(
        sns.cubehelix_palette(
            start=-0.2,
            rot=-0.2,
            n_colors=len(modelling_scenario.pv_module.pv_cells_and_cell_strings),
        )
    )

    # Fix nan errors:
    initial_time = datetime.strptime(
        weather_data[location.name][TIME][0], "%Y-%m-%d %H:%M"
    )
    irradiance_frame = irradiance_frame.fillna(0)
    os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)

    # Determine what type of operation to carry out.
    try:
        operating_mode = OperatingMode(parsed_args.operating_mode)
    except ValueError:
        raise ValueError(
            f"Invalid operating mode specified: {parsed_args.operating_mode}"
        ) from None

    match operating_mode.value:
        case OperatingMode.HOURLY_MPP.value:
            ################################
            # Parallelised MPP computation #
            ################################

            # _, mpp_power, bypassed_cell_strings = process_single_mpp_calculation_without_pbar(
            #     8,
            #     irradiance_frame=irradiance_frame,
            #     locations_to_weather_and_solar_map=locations_to_weather_and_solar_map,
            #     pv_system=pv_system,
            #     scenario=modelling_scenario,
            # )

            # Use joblib to parallelize the for loop
            start_time = time.time()
            print(
                (this_string := "Parallel MPP computation")
                + "." * (88 - (len(this_string) + len(DONE)))
                + " ",
                end="",
            )
            with tqdm(
                desc="MPP computation", total=parsed_args.iteration_length, unit="hour"
            ) as pbar:
                # with ThreadPool(8) as mpool:
                #     results_map = mpool.map(
                #         functools.partial(
                #             process_single_mpp_calculation,
                #             irradiance_frame=irradiance_frame,
                #             locations_to_weather_and_solar_map=locations_to_weather_and_solar_map,
                #             pbar=pbar,
                #             pv_system=pv_system,
                #             scenario=modelling_scenario,
                #         ),
                #         range(parsed_args.start_day_index, parsed_args.start_day_index + parsed_args.iteration_length),
                #     )

                results = Parallel(n_jobs=8)(
                    delayed(
                        functools.partial(
                            process_single_mpp_calculation_without_pbar,
                            irradiance_frame=irradiance_frame,
                            locations_to_weather_and_solar_map=locations_to_weather_and_solar_map,
                            pv_system=pv_system,
                            scenario=modelling_scenario,
                        )
                    )(time_of_day)
                    for time_of_day in range(
                        parsed_args.start_day_index,
                        parsed_args.start_day_index + parsed_args.iteration_length,
                    )
                )

            end_time = time.time()
            print(DONE)
            print(f"Parallel processing took {end_time - start_time:.2f} seconds")

            # Process results and accumulate daily data
            all_mpp_data: list[
                tuple[
                    str,
                    int,
                    dict[BypassedCellString | PVCell, bool],
                    dict[BypassedCellString | PVCell, float],
                ]
            ] = []

            daily_data = defaultdict(list)
            for result in results:
                if result is not None:
                    hour, mpp_power, bypassed_cell_strings, cellwise_mpp = result
                    if mpp_power is not None:
                        daily_data[
                            (
                                date_str := (
                                    initial_time + timedelta(hours=hour)
                                ).strftime("%d/%m/%Y")
                            )
                        ].append((hour, mpp_power))
                        all_mpp_data.append(
                            (
                                date_str,
                                hour,
                                mpp_power,
                                bypassed_cell_strings,
                                cellwise_mpp,
                            )
                        )

            all_mpp_data = [
                [
                    entry[0],
                    entry[1],
                    entry[2],
                    {key.cell_id: bool(value) for key, value in entry[3].items()},
                    {key.cell_id: value for key, value in entry[4].items()},
                ]
                for entry in all_mpp_data
            ]

            import pdb

            pdb.set_trace()

            # Save the output data
            with open(
                os.path.join(
                    OUTPUT_DIRECTORY,
                    f"hourly_mpp_{scenario.name}_"
                    f"{(start_hour:=parsed_args.start_day_index)}_to_"
                    f"{(end_hour:=start_hour + parsed_args.iteration_length)}.json",
                ),
                "w",
                encoding="UTF-8",
            ) as output_file:
                json.dump(all_mpp_data, output_file, indent=4)

            # Generate a heatmap of whether the power generation in the cells.
            cellwise_mpp_frame = pd.DataFrame(
                {entry[1]: entry[4].values() for entry in all_mpp_data}
            ).transpose()
            (cmap_greens := plt.get_cmap("Oranges").copy()).set_under("none")
            (cmap_purples := plt.get_cmap("Blues").copy()).set_under("none")

            plt.figure(figsize=(40 / 5, 32 / 5))
            sns.heatmap(
                cellwise_mpp_frame,
                cmap=cmap_greens,
                cbar_kws={"label": "Power generation in cells / W"},
                square=True,
                vmin=0,
            )
            sns.heatmap(
                -cellwise_mpp_frame,
                cmap=cmap_purples,
                cbar_kws={"label": "Power dissipation in cells / W"},
                square=True,
                vmin=0,
            )
            plt.savefig(
                f"cellwise_mpp_{modelling_scenario.name}_{start_hour}_{end_hour}."
                f"{(fig_format:='pdf')}",
                format=fig_format,
                bbox_inches="tight",
                pad_inches=0,
            )
            # plt.show()

            # Generate a heatmap of whether the cells were in reverse bias.
            reverse_bias_frame = (
                pd.DataFrame({entry[1]: entry[3].values() for entry in all_mpp_data})
                .transpose()
                .astype(bool)
            )
            reverse_bias_frame = (
                reverse_bias_frame.astype(str)
                .replace("0", None)
                .replace("False", 0)
                .replace("false", 0)
                .replace("True", 1)
            )

            sns.heatmap(reverse_bias_frame, cmap="PiYG_r", square=True, vmin=0, vmax=1)
            plt.savefig(
                f"cellwise_reverse_bias_{modelling_scenario.name}_{start_hour}_"
                f"{end_hour}.{(fig_format:='pdf')}",
                format=fig_format,
                bbox_inches="tight",
                pad_inches=0,
            )
            # plt.show()

        case OperatingMode.HOURLY_MPP_MACHINE_LEARNING_TRAING.value:
            #################################
            # BSc. Machine Learning Project #
            #################################

            # Open the data from the training data-set file if it already exists,
            # otherwise, construct new data.
            if os.path.isfile((filename := "training_data_hours.csv")):
                training_data_hours = pd.read_csv(filename)
                # pd.DataFrame
                # 0, 0
                # 1, 0
                # ...
                # 8, 1
                # 9, 0
            else:
                # Compute a series of random hours throughout the year
                # Save them to the file
                pass

            # Split the weather data into training and test data.

            # Use joblib to parallelize the for loop
            start_time = time.time()
            print(
                (this_string := "Parallel MPP computation")
                + "." * (88 - (len(this_string) + len(DONE)))
                + " ",
                end="",
            )
            with tqdm(
                desc="MPP computation", total=parsed_args.iteration_length, unit="hour"
            ) as pbar:
                # with ThreadPool(8) as mpool:
                #     results_map = mpool.map(
                #         functools.partial(
                #             process_single_mpp_calculation,
                #             irradiance_frame=irradiance_frame,
                #             locations_to_weather_and_solar_map=locations_to_weather_and_solar_map,
                #             pbar=pbar,
                #             pv_system=pv_system,
                #             scenario=modelling_scenario,
                #         ),
                #         range(parsed_args.start_day_index, parsed_args.start_day_index + parsed_args.iteration_length),
                #     )

                results = Parallel(n_jobs=8)(
                    delayed(
                        functools.partial(
                            process_single_mpp_calculation_without_pbar,
                            irradiance_frame=irradiance_frame,
                            locations_to_weather_and_solar_map=locations_to_weather_and_solar_map,
                            pv_system=pv_system,
                            scenario=modelling_scenario,
                        )
                    )(time_of_day)
                    for time_of_day in range(
                        parsed_args.start_day_index,
                        parsed_args.start_day_index + parsed_args.iteration_length,
                    )
                )

            end_time = time.time()
            print(DONE)
            print(f"Parallel processing took {end_time - start_time:.2f} seconds")

            # Process results and accumulate daily data
            all_mpp_data: list[
                tuple[
                    str,
                    int,
                    dict[BypassedCellString | PVCell, bool],
                    dict[BypassedCellString | PVCell, float],
                ]
            ] = []

            daily_data = defaultdict(list)
            for result in results:
                if result is not None:
                    hour, mpp_power, bypassed_cell_strings, cellwise_mpp = result
                    if mpp_power is not None:
                        daily_data[
                            (
                                date_str := (
                                    initial_time + timedelta(hours=hour)
                                ).strftime("%d/%m/%Y")
                            )
                        ].append((hour, mpp_power))
                        all_mpp_data.append(
                            (
                                date_str,
                                hour,
                                mpp_power,
                                bypassed_cell_strings,
                                cellwise_mpp,
                            )
                        )

            def _process_bypassing(entry_to_process: bool | None) -> None:
                """
                Process the bypass-diode.

                :param: entry_to_process
                    The entry to process.

                :returns: The processed entry.

                """

                if entry_to_process == "0":
                    return None
                if entry_to_process == "False":
                    return 0
                return 1

            all_mpp_data = [
                [
                    entry[0],
                    entry[1],
                    entry[2],
                    {
                        key.cell_id: _process_bypassing(str(value))
                        for key, value in entry[3].items()
                    },
                    {key.cell_id: value for key, value in entry[4].items()},
                ]
                for entry in all_mpp_data
            ]

            # Save the output data
            with open(
                os.path.join(
                    OUTPUT_DIRECTORY,
                    f"hourly_mpp_{scenario.name}_"
                    f"{(start_hour:=parsed_args.start_day_index)}_to_"
                    f"{(end_hour:=start_hour + parsed_args.iteration_length)}.json",
                ),
                "w",
                encoding="UTF-8",
            ) as output_file:
                json.dump(all_mpp_data, output_file)

            training_data = results[training_hours]

            # Week 3 only: Select the machine-learning algorithm from either an argument
            # that we take in on the command-line, in Python, or from some file, or we
            # go through a list of different algorithms.

            # Run the code, in a parallel loop, to develop a machine-learning object
            # which is trained on the training data.

            # Use the test data to assess how well it worked.

            # Save all the objects and the results.

        case OperatingMode.IRRADIANCE_LINEPLOTS.value:
            ################################
            # Irradiance lineplot plotting #
            ################################

            # Line plot of the irradiance throughout the day.
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
            plt.figure(figsize=(48 / 5, 32 / 5))
            dashes = Dashes()

            for time_of_day in tqdm(
                range(
                    (start_hour := parsed_args.start_day_index),
                    (
                        end_hour := parsed_args.start_day_index
                        + parsed_args.iteration_length
                    ),
                ),
                desc="Plotting irradiance series",
            ):
                sns.lineplot(
                    x=[
                        pv_cell.cell_id
                        for pv_cell in modelling_scenario.pv_module.pv_cells
                    ],
                    y=[
                        1000
                        * irradiance_frame.set_index("hour")
                        .iloc[time_of_day]
                        .values[pv_cell.cell_id]
                        for pv_cell in modelling_scenario.pv_module.pv_cells
                    ],
                    # color="orange",
                    # marker="h",
                    # s=100,
                    # alpha=0.35,
                    dashes=next(dashes),
                    label=_sanitise_time(time_of_day, "%H:%M", initial_time),
                )

            plt.xlabel("Cell ID")
            plt.ylabel("Irradiance / W/m$^2$")
            # plt.legend().remove()
            plt.legend()

            # norm = plt.Normalize(
            #     int(_sanitise_time(start_hour, "%H", initial_time)),
            #     int(_sanitise_time(end_hour - 1, "%H", initial_time)) + 1,
            # )
            # scalar_mappable = plt.cm.ScalarMappable(
            #     cmap=mcolors.LinearSegmentedColormap.from_list(
            #         "Custom", sns.color_palette().as_hex(), parsed_args.iteration_length
            #     ),
            #     norm=norm,
            # )

            # colorbar = (axis := plt.gca()).figure.colorbar(
            #     scalar_mappable,
            #     ax=axis,
            #     label="Hour of the day",
            #     pad=(_pad := 0.025),
            # )
            # plt.title((initial_time + timedelta(hours=plotting_time_of_day)).strftime("%H:%M on %d/%m/%Y"))
            plt.savefig(
                f"irradiance_graph_{modelling_scenario.name}_{start_hour}_{end_hour}."
                f"{(fig_format:='pdf')}",
                format=fig_format,
                bbox_inches="tight",
                pad_inches=0,
            )
            plt.show()

            plt.figure(figsize=(48 / 5, 32 / 5))
            dashes = Dashes()

            for time_of_day in tqdm(
                range(start_hour, end_hour),
                desc="Plotting temperature series",
            ):
                plt.plot(
                    [
                        pv_cell.cell_id
                        for pv_cell in modelling_scenario.pv_module.pv_cells
                    ],
                    [
                        pv_cell.average_cell_temperature(
                            locations_to_weather_and_solar_map[
                                modelling_scenario.location
                            ][time_of_day][TEMPERATURE]
                            + ZERO_CELSIUS_OFFSET,
                            1000
                            * irradiance_frame.set_index("hour")
                            .iloc[time_of_day]
                            .iloc[pv_cell.cell_id],
                            0,
                        )
                        - ZERO_CELSIUS_OFFSET
                        for pv_cell in modelling_scenario.pv_module.pv_cells
                    ],
                    # color="orange",
                    # marker="h",
                    # s=100,
                    # alpha=0.35,
                    dashes=next(dashes),
                    label=_sanitise_time(time_of_day, "%H:%M", initial_time),
                )

            plt.xlabel("Cell ID")
            plt.ylabel("Temperature / Degrees Celsius")
            plt.legend()
            # plt.legend().remove()

            # norm = plt.Normalize(
            #     int(_sanitise_time(start_hour, "%H", initial_time)),
            #     int(_sanitise_time(end_hour - 1, "%H", initial_time)) + 1,
            # )
            # scalar_mappable = plt.cm.ScalarMappable(
            #     cmap=mcolors.LinearSegmentedColormap.from_list(
            #         "Custom", sns.color_palette().as_hex(), parsed_args.iteration_length
            #     ),
            #     norm=norm,
            # )

            # (axis := plt.gca()).figure.colorbar(
            #     scalar_mappable,
            #     ax=axis,
            #     label="Hour of the day",
            #     pad=(_pad := 0.025),
            # )
            # plt.title((initial_time + timedelta(hours=plotting_time_of_day)).strftime("%H:%M on %d/%m/%Y"))
            plt.savefig(
                f"temperature_graph_{modelling_scenario.name}_{start_hour}_{end_hour}."
                f"{(fig_format:='pdf')}",
                format=fig_format,
                bbox_inches="tight",
                pad_inches=0,
            )
            plt.show()

        case OperatingMode.IRRADIANCE_HEATMAPS.value:
            ###############################
            # Irradiance heatmap plotting #
            ###############################

            try:
                scenario_index: int = [
                    index
                    for index, scenario in enumerate(scenarios)
                    if scenario.name == parsed_args.scenario
                ][0]
            except IndexError:
                raise IndexError(
                    "Unable to find current scenario weather information."
                ) from None

            frame_to_plot = cellwise_irradiance_frames[0][1]
            plot_irradiance_with_marginal_means(
                frame_to_plot,
                start_index=(start_hour := parsed_args.start_day_index),
                figname=f"{scenario.name}_{parsed_args.start_day_index}_irradiance",
                heatmap_vmax=(
                    1.0
                    # frame_to_plot.set_index("hour")
                    # .iloc[start_hour : start_hour + 24]
                    # .max()
                    # .max()
                ),
                initial_date=initial_time,
                irradiance_bar_vmax=(
                    1.0
                    # max(
                    #     frame_to_plot.set_index("hour")
                    #     .iloc[start_hour : start_hour + 24]
                    #     .mean(axis=1)
                    #     .max(),
                    #     0,
                    # )
                ),
                irradiance_scatter_vmin=0,
                irradiance_scatter_vmax=(
                    1.25
                    # * (
                    #     max(
                    #         frame_to_plot.set_index("hour")
                    #         .iloc[start_hour : start_hour + 24]
                    #         .mean(axis=0)
                    #         .max(),
                    #         0,
                    #     )
                    # )
                ),
                show_figure=True,
            )

            plot_temperature_with_marginal_means(
                frame_to_plot,
                start_index=(start_hour := parsed_args.start_day_index),
                figname=f"{scenario.name}_{parsed_args.start_day_index}_temperature",
                heatmap_vmax=80,
                initial_date=initial_time,
                locations_to_weather_and_solar_map=locations_to_weather_and_solar_map,
                scenario=scenario,
                show_figure=True,
                temperature_bar_vmax=80,
                temperature_scatter_vmax=80,
            )

        case OperatingMode.IV_CURVES.value:
            ##############################################
            # Plot IV curves for the cells in the module #
            ##############################################

            time_of_day = parsed_args.start_day_index

            current_series = np.linspace(
                0,
                1.1
                * max(
                    [
                        pv_cell.short_circuit_current
                        for pv_cell in scenario.pv_module.pv_cells
                    ]
                ),
                VOLTAGE_RESOLUTION,
            )

            # Maps for plotting curves
            cell_to_current_map: dict[BypassedCellString | PVCell, np.ndarray] = {}
            cell_to_voltage_map: dict[BypassedCellString | PVCell, np.ndarray] = {}
            cell_to_bypass_current_map: dict[
                BypassedCellString | PVCell, np.ndarray
            ] = {}
            cell_to_bypass_voltage_map: dict[
                BypassedCellString | PVCell, np.ndarray
            ] = {}

            # Maps for computing the MPP
            mpp_cell_to_current_map: dict[BypassedCellString | PVCell, np.ndarray] = {}
            mpp_cell_to_voltage_map: dict[BypassedCellString | PVCell, np.ndarray] = {}
            mpp_cell_to_power_map: dict[BypassedCellString | PVCell, np.ndarray] = {}
            mpp_cell_to_bypass_operating: dict[
                BypassedCellString | PVCell, np.ndarray
            ] = {}

            for pv_cell in tqdm(
                modelling_scenario.pv_module.pv_cells_and_cell_strings,
                desc="IV calculation",
                leave=False,
            ):
                (
                    current_series,
                    power_series,
                    voltage_series,
                    bypass_diode_operating,
                ) = pv_cell.calculate_iv_curve(
                    (
                        weather_map := locations_to_weather_and_solar_map[
                            modelling_scenario.location
                        ][time_of_day]
                    )[TEMPERATURE],
                    1000
                    * irradiance_frame.set_index("hour")
                    .iloc[time_of_day]
                    .reset_index(drop=True),
                    weather_map[WIND_SPEED],
                    current_series=current_series,
                )

                mpp_cell_to_current_map[pv_cell] = current_series
                mpp_cell_to_voltage_map[pv_cell] = voltage_series
                mpp_cell_to_power_map[pv_cell] = power_series
                mpp_cell_to_bypass_operating[pv_cell] = bypass_diode_operating

            time_string = (initial_time + timedelta(hours=time_of_day)).strftime(
                "%d_%m_%y_at_%H_%M"
            )

            # Plot the IV curves across the module, along with the power that each
            # bypassed group produces
            plt.figure(figsize=(48 / 5, 32 / 5))
            left_axis = plt.gca()
            right_axis = left_axis.twinx()

            for index, pv_cell in enumerate(
                modelling_scenario.pv_module.pv_cells_and_cell_strings
            ):
                left_axis.plot(
                    mpp_cell_to_voltage_map[pv_cell],
                    mpp_cell_to_current_map[pv_cell],
                    label=pv_cell.cell_id,
                    color=f"C{index}",
                )
                right_axis.plot(
                    mpp_cell_to_voltage_map[pv_cell],
                    mpp_cell_to_power_map[pv_cell],
                    "--",
                    label=pv_cell.cell_id,
                    color=f"C{index}",
                )

            left_axis.set_ylim(bottom=-2.5, top=10)
            right_axis.set_ylim(bottom=-10, top=40)
            left_axis.legend(title="Initial index of cell in string", ncol=4)
            plt.xlabel("Cell-wise, or cell-string-wise, voltage / V")
            left_axis.set_ylabel("Module current / A")
            right_axis.set_ylabel("Cell-wise, or cell-string-wise, power / W")

            plt.savefig(
                f"iv_curves_with_power_{time_string}.pdf",
                format="pdf",
                bbox_inches="tight",
                pad_inches=0,
            )
            plt.show()

            # Plot the IV curves across the module as current and power only.
            plt.figure(figsize=(48 / 5, 32 / 5))

            for index, pv_cell in enumerate(
                modelling_scenario.pv_module.pv_cells_and_cell_strings
            ):
                plt.plot(
                    mpp_cell_to_voltage_map[pv_cell],
                    mpp_cell_to_current_map[pv_cell],
                    label=pv_cell.cell_id,
                    color=f"C{index}",
                )

            plt.ylim(bottom=0, top=10)
            plt.legend(title="Initial index of cell in string", ncol=4)
            plt.xlabel("Cell-wise, or cell-string-wise, voltage / V")
            plt.ylabel("Module current / A")

            plt.savefig(
                f"iv_curves_{time_string}.pdf",
                format="pdf",
                bbox_inches="tight",
                pad_inches=0,
            )
            plt.show()

            # Compute the total power produced
            cell_voltage_interpreters = {
                pv_cell: lambda i: np.interp(
                    i,
                    list(reversed(mpp_cell_to_current_map[pv_cell])),
                    list(reversed(mpp_cell_to_voltage_map[pv_cell])),
                    period=np.max(mpp_cell_to_current_map[pv_cell])
                    - np.min(mpp_cell_to_current_map[pv_cell]),
                )
                for pv_cell in modelling_scenario.pv_module.pv_cells_and_cell_strings
            }

            cellwise_voltage = {
                pv_cell: [
                    cell_voltage_interpreters[pv_cell](current)
                    for current in tqdm(
                        current_series[::CURRENT_SAMPLING_RATE],
                        desc="Interpolating current",
                        leave=False,
                    )
                ]
                for pv_cell in tqdm(
                    modelling_scenario.pv_module.pv_cells_and_cell_strings,
                    desc="String-wise calculation",
                    leave=False,
                )
            }
            module_voltage = [
                sum(sublist) for sublist in zip(*cellwise_voltage.values())
            ]

            for index, pv_cell in enumerate(
                modelling_scenario.pv_module.pv_cells_and_cell_strings
            ):
                plt.plot(
                    mpp_cell_to_voltage_map[pv_cell],
                    mpp_cell_to_current_map[pv_cell],
                    label=pv_cell.cell_id,
                    color=f"C{index}",
                )

            plt.plot(
                module_voltage,
                current_series[::CURRENT_SAMPLING_RATE],
                "--",
                color="orange",
            )

            plt.ylim(bottom=0, top=10)
            plt.legend(title="Initial index of cell in string", ncol=4)
            plt.xlabel("Cell-wise, or cell-string-wise, voltage / V")
            plt.ylabel("Module current / A")

            plt.savefig(
                f"iv_curves_with_module_{time_string}.pdf",
                format="pdf",
                bbox_inches="tight",
                pad_inches=0,
            )
            plt.show()

            module_power = module_voltage * current_series[::CURRENT_SAMPLING_RATE]
            mpp_index: int = list(module_power).index(np.max(module_power))

            plt.figure(figsize=(48 / 5, 32 / 5))
            left_axis = plt.gca()
            right_axis = left_axis.twinx()

            for index, pv_cell in enumerate(
                modelling_scenario.pv_module.pv_cells_and_cell_strings
            ):
                left_axis.plot(
                    mpp_cell_to_current_map[pv_cell],
                    mpp_cell_to_power_map[pv_cell],
                    label=pv_cell.cell_id,
                    color=f"C{index}",
                )

            right_axis.plot(
                current_series[::CURRENT_SAMPLING_RATE],
                module_power,
                "--",
                color="orange",
                label="Module power",
            )
            right_axis.scatter(
                current_series[::CURRENT_SAMPLING_RATE][mpp_index],
                module_power[mpp_index],
                marker="H",
                color="orange",
                s=200,
                label="Maximum power point (MPP)",
            )

            left_axis.set_ylim(bottom=-25, top=25)
            right_axis.set_ylim(
                bottom=-(power_limit := 1.1 * np.max(module_power)), top=power_limit
            )
            plt.xlim(np.min(current_series), np.max(current_series))

            l_handles, l_labels = left_axis.get_legend_handles_labels()
            r_handles, r_labels = right_axis.get_legend_handles_labels()

            left_axis.legend(
                l_handles + [None] + r_handles,
                l_labels + ["Power"] + r_labels,
                #  title="Initial index of cell in string",
                ncol=4,
            )
            plt.xlabel("Module current / A")
            left_axis.set_ylabel("Cell-wise, or cell-string-wise, power / W")
            right_axis.set_ylabel("Module power / W")

            plt.savefig(
                f"ip_curves_{time_string}.pdf",
                format="pdf",
                bbox_inches="tight",
                pad_inches=0,
            )
            plt.show()

            import pdb

            pdb.set_trace()

        case OperatingMode.VALIDATION.value:
            ############################################
            # Validate the model against inputted data #
            ############################################

            # Calcualte the MPP of the PV system for each of the horus to be modelled.
            if parsed_args.timestamps_file is None:
                raise ArgumentError(
                    "Cannot carry out hourly calculation without timestamps file."
                )

            # Open the timestamps file.
            try:
                with open(parsed_args.timestamps_file, "r") as f:
                    timestamps_data = pd.read_csv(f, index_col="timestamp")
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"Could not find timestamps file: {parsed_args.timestamps_file}"
                )

            # Compute the start time within the year based on the time stamp if not provided.
            if (_start_time_column_name := "start_time") not in timestamps_data.columns:
                timestamps_data[_start_time_column_name] = [
                    (
                        datetime.datetime.combine(
                            datetime.date(
                                year=int(entry[0].split(" ")[0].split("/")[2]),
                                month=int(entry[0].split(" ")[0].split("/")[1]),
                                day=int(entry[0].split(" ")[0].split("/")[0]),
                            ),
                            datetime.time(
                                hour=int(entry[0].split(" ")[1].split(":")[0])
                            ),
                        )
                        - datetime.datetime(year=2023, month=1, day=1, hour=0, minute=0)
                    ).total_seconds()
                    / 3600
                    for entry in timestamps_data.iterrows()
                ]

            # For each hour within the series of start times, compute the MPP
            # unless the file for validation already exists with these times in.
            if os.path.isfile(
                validation_filename := "validation_file_{timestamps_filename}_{start_time}_{end_time}".format(
                    timestamps_filename=parsed_args.timestamps_file,
                    start_time=int(timestamps_data[_start_time_column_name][0]),
                    end_time=int(timestamps_data[_start_time_column_name][-1]),
                )
            ):
                with open(
                    validation_filename, "r", encoding="UTF-8"
                ) as validation_file:
                    all_mpp_data = json.load(validation_file)

            else:
                results = Parallel(n_jobs=8)(
                    delayed(
                        functools.partial(
                            process_single_mpp_calculation_without_pbar,
                            irradiance_frame=irradiance_frame,
                            locations_to_weather_and_solar_map=locations_to_weather_and_solar_map,
                            pv_system=pv_system,
                            scenario=modelling_scenario,
                        )
                    )(int(time_of_day))
                    for time_of_day in timestamps_data[_start_time_column_name]
                )

                # Format the data correctly.
                all_mpp_data: list[
                    tuple[
                        str,
                        int,
                        dict[BypassedCellString | PVCell, bool],
                        dict[BypassedCellString | PVCell, float],
                    ]
                ] = []

                daily_data = defaultdict(list)
                for result in results:
                    if result is not None:
                        hour, mpp_power, bypassed_cell_strings, cellwise_mpp = result
                        if mpp_power is not None:
                            daily_data[
                                (
                                    date_str := (
                                        initial_time + timedelta(hours=hour)
                                    ).strftime("%d/%m/%Y")
                                )
                            ].append((hour, mpp_power))
                            all_mpp_data.append(
                                (
                                    date_str,
                                    hour,
                                    mpp_power,
                                    bypassed_cell_strings,
                                    cellwise_mpp,
                                )
                            )

                all_mpp_data = [
                    [
                        entry[0],
                        entry[1],
                        entry[2],
                        {key.cell_id: bool(value) for key, value in entry[3].items()},
                        {key.cell_id: value for key, value in entry[4].items()},
                    ]
                    for entry in all_mpp_data
                ]

                with open(
                    validation_filename, "w", encoding="UTF-8"
                ) as validation_file:
                    json.dumps(all_mpp_data, validation_file, indent=4)

            import pdb

            pdb.set_trace()

        case _:
            raise ArgumentError(
                "The operating mode specified is not implemented. Exiting."
            )

    # # Save daily data to a single Excel file with separate sheets
    # with pd.ExcelWriter(
    #     os.path.join(
    #         OUTPUT_DIRECTORY, f"mpp_daily_summary_{modelling_scenario.name}.xlsx"
    #     ),
    #     engine="openpyxl",
    # ) as writer:

    #     # Ensure at least one sheet is present and visible
    #     workbook = writer.book
    #     placeholder_sheet = workbook.create_sheet(title="Sheet1")

    #     for date_str, data in daily_data.items():
    #         df = pd.DataFrame(data, columns=[HOUR, "Power / W"])
    #         total_mpp = df["Power / W"].sum()
    #         df.loc["Total"] = ["Total", total_mpp]  # Adding the total MPP at the end
    #         df.to_excel(writer, sheet_name=date_str, index=False)

    #     # Remove the placeholder sheet if it was not used
    #     if "Sheet1" in workbook.sheetnames and len(workbook.sheetnames) > 1:
    #         del workbook["Sheet1"]

    # # #TODO:
    # # # - Improve the speed of the calculation so it can be run for all hours.
    # # # - Some way to store whether the cells have been bypassed.

    # def _post_process_split_axes(ax1, ax2):
    #     """
    #     Function to post-process the joining of axes.
    #     Adapted from:
    #         https://matplotlib.org/stable/gallery/subplots_axes_and_figures/broken_axis.html
    #     """
    #     # hide the spines between ax and ax2
    #     ax1.spines.bottom.set_visible(False)
    #     ax1.spines.top.set_visible(False)
    #     ax2.spines.top.set_visible(False)
    #     ax1.tick_params(
    #         labeltop=False, labelbottom=False
    #     )  # don't put tick labels at the top
    #     ax2.xaxis.tick_bottom()
    #     # Now, let's turn towards the cut-out slanted lines.
    #     # We create line objects in axes coordinates, in which (0,0), (0,1),
    #     # (1,0), and (1,1) are the four corners of the axes.
    #     # The slanted lines themselves are markers at those locations, such that the
    #     # lines keep their angle and position, independent of the axes size or scale
    #     # Finally, we need to disable clipping.
    #     d = 0.5  # proportion of vertical to horizontal extent of the slanted line
    #     kwargs = dict(
    #         marker=[(-1, -d), (1, d)],
    #         markersize=12,
    #         linestyle="none",
    #         color="k",
    #         mec="k",
    #         mew=1,
    #         clip_on=False,
    #     )
    #     ax1.plot([0], [0], transform=ax1.transAxes, **kwargs)
    #     ax2.plot([0], [1], transform=ax2.transAxes, **kwargs)

    # gridspec = {"hspace": 0.1, "height_ratios": [1, 1, 0.4, 1, 1]}
    # fig, axes = plt.subplots(5, 2, figsize=(48 / 5, 32 / 5), gridspec_kw=gridspec)
    # fig.subplots_adjust(hspace=0, wspace=0.25)

    # axes[2, 0].set_visible(False)
    # axes[2, 1].set_visible(False)
    # y_label_coord: int = int(-850)

    # axes[0, 0].get_shared_x_axes().join(axes[0, 0], axes[1, 0])
    # axes[3, 0].get_shared_x_axes().join(axes[3, 0], axes[4, 0])
    # axes[3, 1].get_shared_x_axes().join(axes[3, 1], axes[4, 1])
    # axes[0, 1].get_shared_x_axes().join(axes[0, 1], axes[1, 1])

    # curve_info = pvlib.pvsystem.singlediode(
    #     photocurrent=IL,
    #     saturation_current=I0,
    #     resistance_series=Rs,
    #     resistance_shunt=Rsh,
    #     nNsVth=nNsVth,
    #     ivcurve_pnts=100,
    #     method="lambertw",
    # )
    # plt.plot(curve_info["v"], curve_info["i"])
    # plt.show()

    # pvlib.singlediode.bishop88_i_from_v(
    #     -14.95,
    #     photocurrent=IL,
    #     saturation_current=I0,
    #     resistance_series=Rs,
    #     resistance_shunt=Rsh,
    #     nNsVth=nNsVth,
    #     breakdown_voltage=-15,
    #     breakdown_factor=2e-3,
    #     breakdown_exp=3,
    # )

    # v_oc = pvlib.singlediode.bishop88_v_from_i(
    #     0.0,
    #     photocurrent=IL,
    #     saturation_current=I0,
    #     resistance_series=Rs,
    #     resistance_shunt=Rsh,
    #     nNsVth=nNsVth,
    #     method="lambertw",
    # )
    # voltage_array = np.linspace(-15 * 0.999, v_oc, 1000)
    # ivcurve_i, ivcurve_v, _ = pvlib.singlediode.bishop88(
    #     voltage_array,
    #     photocurrent=IL,
    #     saturation_current=I0,
    #     resistance_series=Rs,
    #     resistance_shunt=Rsh,
    #     nNsVth=nNsVth,
    #     breakdown_voltage=-15,
    # )

    # # Determine the scenario index

    # frame_slice = (
    #     cellwise_irradiance_frames[scenario_index][1]
    #     .iloc[(start_index := parsed_args.start_day_index) : start_index + 24]
    #     .set_index("hour")
    # )
    # sns.heatmap(
    #     frame_slice,
    #     cmap=sns.blend_palette(
    #         [
    #             "#144E56",
    #             "#28769C",
    #             "teal",
    #             "#94B49F",
    #             "grey",
    #             "silver",
    #             "orange",
    #             "#E04606",
    #         ],
    #         as_cmap=True,
    #     ),
    #     vmin=0,
    #     cbar_kws={"label": "Irradiance / kWm$^{-2}$"},
    # )
    # plt.xlabel("Cell index within panel")
    # plt.ylabel("Hour of the day")
    # plt.show()

    # import pdb

    # pdb.set_trace()

    # plot_irradiance_with_marginal_means(
    #     cellwise_irradiance_frames[1][1],
    #     start_index=start_index,
    #     figname="july_eighth_medium_panel",
    #     heatmap_vmax=heatmap_vmax,
    #     irradiance_bar_vmax=irradiance_bar_vmax,
    #     irradiance_scatter_vmin=0,
    #     irradiance_scatter_vmax=irradiance_scatter_vmax,
    # )
    # plot_irradiance_with_marginal_means(
    #     cellwise_irradiance_frames[2][1],
    #     start_index=start_index,
    #     figname="july_eighth_large_panel",
    #     heatmap_vmax=heatmap_vmax,
    #     irradiance_bar_vmax=irradiance_bar_vmax,
    #     irradiance_scatter_vmin=0,
    #     irradiance_scatter_vmax=irradiance_scatter_vmax,
    # )

    # plot_irradiance_with_marginal_means(
    #     cellwise_irradiance_frames[0][1],
    #     start_index=(start_index := 24 * 31 * 7 + 48),
    #     figname="august_eight_small_panel",
    #     heatmap_vmax=(
    #         heatmap_vmax := cellwise_irradiance_frames[2][1]
    #         .set_index("hour")
    #         .iloc[start_index : start_index + 24]
    #         .max()
    #         .max()
    #     ),
    #     irradiance_bar_vmax=(
    #         irradiance_bar_vmax := (
    #             max(
    #                 cellwise_irradiance_frames[0][1]
    #                 .set_index("hour")
    #                 .iloc[start_index : start_index + 24]
    #                 .mean(axis=1)
    #                 .max(),
    #                 cellwise_irradiance_frames[1][1]
    #                 .set_index("hour")
    #                 .iloc[start_index : start_index + 24]
    #                 .mean(axis=1)
    #                 .max(),
    #                 cellwise_irradiance_frames[2][1]
    #                 .set_index("hour")
    #                 .iloc[start_index : start_index + 24]
    #                 .mean(axis=1)
    #                 .max(),
    #             )
    #         )
    #     ),
    #     irradiance_scatter_vmin=0,
    #     irradiance_scatter_vmax=(
    #         irradiance_scatter_vmax := 1.25
    #         * (
    #             max(
    #                 cellwise_irradiance_frames[0][1]
    #                 .set_index("hour")
    #                 .iloc[start_index : start_index + 24]
    #                 .sum(axis=0)
    #                 .max(),
    #                 cellwise_irradiance_frames[1][1]
    #                 .set_index("hour")
    #                 .iloc[start_index : start_index + 24]
    #                 .sum(axis=0)
    #                 .max(),
    #                 cellwise_irradiance_frames[2][1]
    #                 .set_index("hour")
    #                 .iloc[start_index : start_index + 24]
    #                 .sum(axis=0)
    #                 .max(),
    #             )
    #         )
    #     ),
    # )
    # plot_irradiance_with_marginal_means(
    #     cellwise_irradiance_frames[1][1],
    #     start_index=start_index,
    #     figname="august_eighth_medium_panel",
    #     heatmap_vmax=heatmap_vmax,
    #     irradiance_bar_vmax=irradiance_bar_vmax,
    #     irradiance_scatter_vmin=0,
    #     irradiance_scatter_vmax=irradiance_scatter_vmax,
    # )
    # plot_irradiance_with_marginal_means(
    #     cellwise_irradiance_frames[2][1],
    #     start_index=start_index,
    #     figname="august_eighth_large_panel",
    #     heatmap_vmax=heatmap_vmax,
    #     irradiance_bar_vmax=irradiance_bar_vmax,
    #     irradiance_scatter_vmin=0,
    #     irradiance_scatter_vmax=irradiance_scatter_vmax,
    # )

    # weather = locations_with_weather_and_solar[scenario.location]
    # data_to_scatter = weather.iloc[start_index : start_index + 24]
    # plt.scatter(
    #     range(24),
    #     data_to_scatter["irradiance_direct"],
    #     marker="H",
    #     s=40,
    #     label="Direct irradiance",
    # )
    # plt.scatter(
    #     range(24),
    #     data_to_scatter["irradiance_diffuse"],
    #     marker="H",
    #     s=40,
    #     label="Diffuse irradiance",
    # )
    # plt.scatter(
    #     range(24),
    #     data_to_scatter["irradiance_direct_normal"],
    #     marker="H",
    #     s=40,
    #     label="DNI irradiance",
    # )
    # plt.legend()
    # plt.show()

    # (
    #     frame_slice := cellwise_irradiance_frames[2][1].iloc[
    #         start_index : start_index + 24
    #     ]
    # ).plot(
    #     x="hour",
    #     y=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
    #     colormap=sns.color_palette("PuBu_r", n_colors=18, as_cmap=True),
    # )


if __name__ == "__main__":
    main(sys.argv[1:])
