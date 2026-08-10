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

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from rich.progress import track

from src.polytunnelpv.__utils__ import (
    BOLTZMAN_CONSTANT,
    ELECTRON_CHARGE,
    VOLTAGE_RESOLUTION,
)
from src.polytunnelpv.pv_module.pv_cell import PVCell, ZERO_CELSIUS_OFFSET

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

    bypass_voltage: float | None
    end_index: int
    start_index: int
    saturation_current: float = 1e-12
    ideality_factor: float = 1.1
    thermal_voltage: float = field(init=False)

    def __repr__(self) -> str:
        """
        Return a nice-looking representation of the bypass diode.

        """

        return (
            f"BypassDiode(bypass_voltage={self.bypass_voltage}, "
            f"start_index={self.start_index}), end_index={self.end_index}, "
            f"saturation_current={self.saturation_current}, "
            f"ideality_factor={self.ideality_factor})"
        )

    def calculate_i_from_v(self, voltage: float | list[float]) -> np.array:
        """
        Calculate the current based on the voltage using the Shockley diode equation.

        :param: voltage
            Either a single value or a list of values for which to calculate.

        :returns: The current.

        """

        voltage = np.array(voltage)  # Ensure voltage is an array for vector operations
        current = self.saturation_current * (
            np.exp(-voltage / (self.ideality_factor * self.thermal_voltage)) - 1
        )
        return current

    def calculate_v_from_i(self, current: float | list[float]) -> np.array:
        """
        Calculate the voltage based on the current using the Shockley diode equation.

        :param: current
            Either a single value or a list of values for which to calculate.

        :returns: The voltage.

        """

        current = np.array(current)
        voltage = -(self.ideality_factor * self.thermal_voltage) * np.log(
            1 + current / self.saturation_current
        )

        return voltage


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

    @property
    def breakdown_voltage(self) -> float:
        """
        Return the breakdown voltage, i.e., the bypass-diode voltage.

        Returns:
            - The breakdown voltage.

        """

        return self.bypass_diode.bypass_voltage

    @property
    def cell_id(self) -> float:
        """
        Return the cell ID for the fist cell in the string.

        Returns:
            - The cell ID of the first cell within the string.

        """

        return min([cell.cell_id for cell in self.pv_cells])

    def __hash__(self) -> int:
        """
        Return a hash of the first cell ID.

        """

        return hash(self.cell_id)

    def set_bypass_diode_temperature(
        self, ambient_celsius_temperature: float, irradiance_array: float
    ) -> None:
        """
        Initialize the bypass diode parameters.

        :param: ambient_celsius_temperature
            The ambient temperature in degC.

        :param: irradiance_array
            The irradiance falling on the system.

        """

        # Calculate the average cell temperature for the set of cells being bypassed.
        temperatures = [
            pv_cell.average_cell_temperature(
                ambient_celsius_temperature + ZERO_CELSIUS_OFFSET,
                irradiance_array[pv_cell.cell_id],
                0,  # Assuming wind speed of 0 for simplicity
            )
            for pv_cell in self.pv_cells
        ]
        average_temperature = np.mean(temperatures)

        # Calculate thermal voltage (kT/q) in Volts
        thermal_voltage = BOLTZMAN_CONSTANT * average_temperature / ELECTRON_CHARGE
        self.bypass_diode.thermal_voltage = thermal_voltage

    def calculate_iv_curve(
        self,
        ambient_celsius_temperature: float,
        irradiance_array: np.ndarray,
        wind_speed: float,
        *,
        current_density_series: np.ndarray | None = None,
        current_series: np.ndarray | None = None,
        voltage_series: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[bool]]:
        """
        Calculate the IV curve for the bypassed string of cells.

        :param: ambient_celsius_temperature
            The ambient temperature, in degrees Celsius.

        :param: irradiance_array
            The irradiance, in W/m^2, striking all the cells in the module.

        :param: wind_speed
            The wind speed, in m/s, over the cell.

        :param: current_density_series
            If provided, the current-density series---the series of opints over which to
            calculate the current an power output from the cell.

        :param: current_series
            The series of current points over which to calculate the current and power
            output from the cell.

        :param: voltage_series
            The series of voltage points over which to calculate the current and power
            output from the cell.

        :returns: current_series
            The current values.

        :returns: power_series
            The power values.

        :returns: voltage_series
            The voltage series.

        :returns: pv_cells_bypassed
            Whether the PV cells have been bypassed.

        """

        # Set the bypass-diode tepmerature.
        self.set_bypass_diode_temperature(ambient_celsius_temperature, irradiance_array)

        # Calculate the curves for each cell
        cell_to_iv_series: dict[PVCell, tuple[np.ndarray, np.ndarray, np.ndarray]] = {
            pv_cell: pv_cell.calculate_iv_curve(
                ambient_celsius_temperature,
                irradiance_array,
                wind_speed,
                current_density_series=current_density_series,
                current_series=current_series,
                # voltage_series=voltage_series,
            )
            for pv_cell in self.pv_cells
            # for pv_cell in track(
            #     self.pv_cells, description="Bypassed IV curves", transient=True
            # )
        }

        # Ensure that the monatomic adjustment to the voltage and current series produces
        # the same current series for each cell and remove .
        monatomic_filter: pd.Series = pd.Series([True] * len(current_series))
        for pv_cell in self.pv_cells:
            monatomic_filter = (
                (_voltage_series := pd.Series(cell_to_iv_series[pv_cell][2]))
                <= pd.Series(_voltage_series).cummin()
            ) & monatomic_filter

        # Apply the monatomic filters.
        for pv_cell in self.pv_cells:
            _filtered_i = cell_to_iv_series[pv_cell][0][monatomic_filter]
            _filtered_p = cell_to_iv_series[pv_cell][1][monatomic_filter]
            _filtered_v = cell_to_iv_series[pv_cell][2][monatomic_filter]
            _filtered_b = np.array(cell_to_iv_series[pv_cell][3])[monatomic_filter]
            cell_to_iv_series[pv_cell] = [
                _filtered_i,
                _filtered_p,
                _filtered_v,
                _filtered_b,
            ]

        # Add up the voltage for each cell
        combined_voltage_series = sum(
            cell_to_iv_series[pv_cell][2] for pv_cell in self.pv_cells
        )

        # Determine the bypass-diode curve
        bypass_diode_curve = self.bypass_diode.calculate_i_from_v(
            bypass_voltage_series := np.linspace(
                min(
                    10 * self.bypass_diode.bypass_voltage, min(combined_voltage_series)
                ),
                max(combined_voltage_series),
                VOLTAGE_RESOLUTION,
            )
        )

        # Interpolate the current series based on the combined voltage series
        cell_current_series = list(cell_to_iv_series.values())[0][0]
        cell_current_curve = np.interp(
            bypass_voltage_series,
            list(reversed(combined_voltage_series)),
            list(reversed(cell_current_series)),
        )

        # Sum up the two current curves
        total_current = cell_current_curve + bypass_diode_curve

        # Compute the power across the module
        total_power = total_current * bypass_voltage_series

        # Determine whether the cells have been bypassed
        pv_cells_bypassed: list[bool] = list(bypass_diode_curve > cell_current_curve)

        # Debug plotting code to visualise IV curves.
        # import matplotlib.pyplot as plt
        # import seaborn as sns

        # plt.figure(figsize=(171 / 25.4, 120 / 25.4))
        # plt.plot(
        #     bypass_voltage_series,
        #     bypass_diode_curve,
        #     color="orange",
        #     label="Bypass diode",
        # )
        # plt.plot(
        #     bypass_voltage_series, cell_current_curve, color="blue", label="PV cell"
        # )
        # plt.plot(bypass_voltage_series, total_current, color="purple", label="Overall")
        # plt.ylim(-0.5, 5)
        # plt.xlabel("Voltage / V")
        # plt.ylabel("Current / A")
        # plt.savefig(f"Bypass diode demonstration: {self.cell_id}.pdf", format="pdf", bbox_inches="tight", pad_inches=0.05)

        return (
            total_current,
            total_power,
            bypass_voltage_series,
            pv_cells_bypassed,
        )

    def calculate_iv_curve_for_plotting(
        self,
        ambient_celsius_temperature: float,
        irradiance_array: np.ndarray,
        wind_speed: float,
        *,
        current_density_series: np.ndarray | None = None,
        current_series: np.ndarray | None = None,
        voltage_series: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[bool]]:
        """
        Calculate the IV curve for the bypassed string of cells for plotting.

        This function differs from the normal calculation in that the bypass-diode curve
        is computed as a function of current within a sensible range, i.e., up to 1.25
        times the short-circuit current.

        :param: ambient_celsius_temperature
            The ambient temperature, in degrees Celsius.

        :param: irradiance_array
            The irradiance, in W/m^2, striking all the cells in the module.

        :param: wind_speed
            The wind speed, in m/s, over the cell.

        :param: current_density_series
            If provided, the current-density series---the series of opints over which to
            calculate the current an power output from the cell.

        :param: current_series
            The series of current points over which to calculate the current and power
            output from the cell.

        :param: voltage_series
            The series of voltage points over which to calculate the current and power
            output from the cell.

        :returns: current_series
            The current values.

        :returns: voltage_series
            The voltage series.

        :returns: bypass_diode_curve
            The curve for the bypass diodes.

        """

        # Set the bypass-diode tepmerature.
        self.set_bypass_diode_temperature(ambient_celsius_temperature, irradiance_array)

        # Calculate the curves for each cell
        cell_to_iv_series: dict[PVCell, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        # for pv_cell in track(
        #     self.pv_cells, description="Bypassed IV curves", transient=True
        # ):
        for pv_cell in self.pv_cells:
            cell_to_iv_series[pv_cell] = pv_cell.calculate_iv_curve(
                ambient_celsius_temperature,
                irradiance_array,
                wind_speed,
                current_density_series=current_density_series,
                current_series=current_series,
                voltage_series=voltage_series,
            )

        # Add up the voltage for each cell
        combined_voltage_series = sum(
            cell_to_iv_series[pv_cell][2] for pv_cell in self.pv_cells
        )

        # Determine the bypass-diode curve
        bypass_diode_curve = self.bypass_diode.calculate_i_from_v(
            (
                bypass_diode_voltage := np.linspace(
                    self.bypass_diode.calculate_v_from_i(np.max(current_series)),
                    np.max(combined_voltage_series),
                    VOLTAGE_RESOLUTION,
                )
            )
        )

        return (
            current_series,
            combined_voltage_series,
            bypass_diode_curve,
            bypass_diode_voltage,
        )
