#!/usr/bin/python3.10
########################################################################################
# pv_system.py - Curved-PV-module system Python module.                                #
#                                                                                      #
# Author: Ben Winchester                                                               #
# Copyright: Ben Winchester, 2024                                                      #
# Date created: 30/05/2024                                                             #
# License: Open source                                                                 #
# Time created: 13:33:00                                                               #
########################################################################################
"""
pv_system.py - The pv-system module for Polytunnel-PV.

This module provides functionality for the modelling of the PV module system: how the
modules (modelled elsewhere) are connected in strings.

"""

from dataclasses import dataclass

import numpy as np

__all__ = (
    "ModuleString",
    "PVSystem",
)


@dataclass
class ModuleString:
    """
    Represents a string of modules.

    .. attribute:: n_modules:
        The number of modules contained within the string.

    """

    n_modules: int


@dataclass
class PVSystem:
    """
    Represents a series of module strings connected into a system.

    .. attributes:: strings
        The module strings in the system.

    """

    strings: list[ModuleString]

    def __post_init__(self) -> None:
        """
        Enforce that the strings are balanced.

        Raises:
            - NotImplementedError:
                Raised if the strings are not balanced.

        """

        if len({string.n_modules for string in self.strings}) > 1:
            raise NotImplementedError(
                "PV systems can only consist of strings of even length."
            )

    def combine_currents(
        self, module_current: float | np.ndarray
    ) -> float | np.ndarray:
        """
        Combine the current across the strings.

        Current sources will add across each string.

        :param: **module_current:**
            The module current characterstic values.

        :returns:
            The combined current for the system.

        """

        return module_current * len(self.strings)

    def combine_voltages(
        self, module_voltage: float | np.ndarray
    ) -> float | np.ndarray:
        """
        Combine voltages across the strings.

        :param: **module_voltage:**
            The module voltage charactersic values.

        :returns:
            The combined voltage for the system.

        """

        # Compue the voltage for each of the strings
        string_voltage: float | np.ndarray = self.strings[0].n_modules * module_voltage

        # Add these voltages together in parallel.
        # NOTE: This will only work if the voltages are equal. Otherwise, an error
        # should be raised.
        if len({string.n_modules for string in self.strings}) > 1:
            raise NotImplementedError(
                "PV systems can only consist of strings of even length."
            )
        return string_voltage

    def combine_powers(self, module_power: float | np.ndarray) -> float | np.ndarray:
        """
        Combine the power across the modules.

        NOTE: Each module will add its own power. So, the total number of modules will
        simply be the total power in the system provided that the configuration is
        valid.

        :param: **module_power:**
            The module power characteristic values.

        :returns:
            The combined power for the system.

        """

        if len({string.n_modules for string in self.strings}) > 1:
            raise NotImplementedError(
                "PV systems can only consist of strings of even length."
            )

        return self.strings[0].n_modules * len(self.strings) * module_power
