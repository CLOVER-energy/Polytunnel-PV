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

__version__ = "1.0.0a1"

import os
import pvlib
import sys
import yaml

from .pv_module.pv_module import (
    CircularCurve,
    Curve,
    CurveType,
    CurvedPVModule,
    ModuleType,
    TYPE_TO_CURVE_MAPPING,
)

# FILE_ENCODING:
#   The encoding to use when opening and closing files.
FILE_ENCODING: str = "UTF-8"

# INPUT_DATA_DIRECTORY:
#   The name of the input-data directory.
INPUT_DATA_DIRECTORY: str = "input_data"

# LOCATIONS_FILENAME:
#   The filename for the locations file.
LOCATIONS_FILENAME: str = "locations.yaml"

# POLYTUNNEL_CURVE@
#   Keyword used for parsing the information about the curve on which a solar panel
# sits.
POLYTUNNEL_CURVE: str = "polytunnel_curve"

# POLYTUNNELS_FILENAME:
#   The name of the polytunnels file.
POLYTUNNELS_FILENAME: str = "polytunnels.yaml"

# PV_MODULES_FILENAME:
#   The name of the PV-modules file.
PV_MODULES_FILENAME: str = "pv_modules.yaml"

# TYPE:
#   Keyword used to determine the module type of the PV.
TYPE: str = "type"


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


def _parse_pv_modules(polytunnels: dict[str, Curve]) -> list[CurvedPVModule]:
    """
    Parse the curved PV module information from the files.

    Inputs:
        - polytunnels:
            A mapping between polytunnel names and instances.

    Outputs:
        The parsed PV modules as a list.

    """

    with open(
        os.path.join(INPUT_DATA_DIRECTORY, PV_MODULES_FILENAME),
        "r",
        encoding=FILE_ENCODING,
    ) as f:
        pv_module_data = yaml.safe_load(f)

    def _construct_pv_module(pv_module_entry) -> CurvedPVModule:
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

        return constructor(**pv_module_entry)

    return [_construct_pv_module(pv_module_entry) for pv_module_entry in pv_module_data]


def main(unparsed_arguments) -> None:
    """
    Main method for Polytunnel-PV.

    """

    locations = _parse_locations()
    polytunnels = _parse_polytunnel_curves()
    pv_modules = _parse_pv_modules(
        {polytunnel.name: polytunnel for polytunnel in polytunnels}
    )

    import pdb

    pdb.set_trace()


if __name__ == "__main__":
    main(sys.argv[1:])
