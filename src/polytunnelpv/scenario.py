#!/usr/bin/python3.10
########################################################################################
# scenario.py - Scenario module for Polytunnel-PV.                                     #
#                                                                                      #
# Author: Ben Winchester                                                               #
# Copyright: Ben Winchester, 2024                                                      #
# Date created: 22/02/2024                                                             #
# License: Open source                                                                 #
########################################################################################
"""
scenario.py - The scenario module for Polytunnel-PV.

"""

from dataclasses import dataclass
from typing import Type, TypeVar

import pvlib

from .pv_module.pv_module import CurvedPVModule
from .__utils__ import NAME

__all__ = "Scenario"

# LOCATION:
#   Keyword for parsing location information.
LOCATION: str = "location"

# PV_MODULE:
#    Keyword for parsing the pv-module information.
PV_MODULE: str = "pv_module"

# Type variable for Curve and children.
S = TypeVar(
    "S",
    bound="Scenario",
)


@dataclass
class Scenario:
    """
    Represents a scenario to run.

    .. attribute:: location
        The location for the scenario.

    .. attribute:: name
        The name of the scenario.

    .. attribute:: pv_module
        The PV Module to use.

    """

    location: pvlib.location.Location
    name: str
    pv_module: CurvedPVModule

    @classmethod
    def from_scenarios_file(
        cls: Type[S],
        entry: dict[str, str],
        locations: dict[str, pvlib.location.Location],
        pv_modules: dict[str, CurvedPVModule],
    ) -> S:
        """
        Instantiate based on data from the scenarios file.

        :param: **entry:**
            The data from the file to use.

        :param: **locations:**
            The `dict` of locations to use, keyed by name.

        :param: **pv_modules:**
            The `dict` of curved PV modules to use, keyed by.

        """

        try:
            return cls(
                locations[entry[LOCATION]], entry[NAME], pv_modules[entry[PV_MODULE]]
            )
        except KeyError:
            raise KeyError("Missing information in scenarios file.") from None

    def __eq__(self, other) -> bool:
        """Says two scenarios are equal if their names are equal."""

        return bool(self.name == other.name)
