#!/usr/bin/python3.10
########################################################################################
# bypass_diode.py - Module to represent bypass diodes.                                 #
#                                                                                      #
# Author: Ben Winchester                                                               #
# Copyright: Ben Winchester, 2024                                                      #
# Date created: 18/05/2024                                                             #
# License: Open source                                                                 #
# Time created: 14:24:00                                                               #
########################################################################################
"""
bypass_diode.py - The bypass-diodde module for Polytunnel-PV.

This module provides functionality for the modelling of bypass diodes within PV modules.

"""

from dataclasses import dataclass

from .pv_cell import PVCell

__all__ = ("BypassDiode", "BypassedCellString")


@dataclass(kw_only=True)
class BypassDiode:
    """
    Represents a bypass diode.

    .. attribute:: bypass_voltage
        The voltage at which the bypass diode will kick in and bypass the cell series.

    .. attribute:: end_index
        The end index for which to bypass cells.

    .. attribute:: start_index
        The start index for which to bypass cells.

    """

    bypass_voltage: float
    end_index: int
    start_index: int


@dataclass(kw_only=True)
class BypassedCellString:
    """
    Represents a series of cells in a string, bypassed by a single diode.

    .. attribute:: bypass_diode
        The bypass diode installed.

    .. attribute:: pv_cells
        The `list` of PV cells that are in series and bypassed by the diode.

    """

    bypass_diode: BypassDiode
    pv_cells: list[PVCell]
