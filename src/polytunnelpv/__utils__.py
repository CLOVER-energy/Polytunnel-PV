#!/usr/bin/python3.10
########################################################################################
# __utils__.py - Utility module for Polytunnel-PV.                                     #
#                                                                                      #
# Author: Ben Winchester                                                               #
# Copyright: Ben Winchester, 2024                                                      #
# Date created: 22/02/2024                                                             #
# License: Open source                                                                 #
########################################################################################
"""
__utils__.py - The utility module for Polytunnel-PV.

"""

__all__ = (
    "BOLTZMAN_CONSTANT",
    "ELECTRON_CHARGE",
    "NAME",
    "VOLTAGE_RESOLUTION",
)

# BOLTZMAN_CONSTANT:
#   The Boltzman constant.
BOLTZMAN_CONSTANT: float = 1.38064852e-23

# ELECTRON_CHARGE:
#    The electron charge, in Coulombs.
ELECTRON_CHARGE: float = 1.60217662e-19

# NAME:
#   Keyword for parsing the name.
NAME: str = "name"

# VOLTAGE_RESOLUTION:
#   The number of points to use for IV plotting.
VOLTAGE_RESOLUTION: int = 1000
