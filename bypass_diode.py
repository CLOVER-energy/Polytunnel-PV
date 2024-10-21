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
bypass_diode.py - The bypass-diode module for Polytunnel-PV.

This module provides functionality for the modeling of bypass diodes within PV modules.

"""

from dataclasses import dataclass, field

import numpy as np

from tqdm import tqdm

from .pv_cell import PVCell, ZERO_CELSIUS_OFFSET

__all__ = ("BypassDiode", "BypassedCellString")


@dataclass(kw_only=True)
class BypassDiode:
    """
    Represents a bypass diode.
    """

    bypass_voltage: float
    end_index: int
    start_index: int
    saturation_current: float = 1e-12
    ideality_factor: float = 1.1
    thermal_voltage: float = field(init=False)

    def calculate_i_from_v(self, voltage: float | list[float]) -> np.array:
        """
        Calculate the current based on the voltage using the Shockley diode equation.
        """
        # import pdb
        # pdb.set_trace()
        voltage = np.array(voltage)  # Ensure voltage is an array for vector operations
        current = self.saturation_current * (
            np.exp(-voltage / (self.ideality_factor * self.thermal_voltage)) - 1
        )
        return current

    # def calculate_i_from_v(self, voltage: float | list[float]):
    #     return 0.5


def calculate_total_current(voltage_series, total_resistance):
    """
    Calculate the total current through the parallel combination of cell and diode.
    """
    voltage_series = np.array(voltage_series)
    total_resistance = np.array(total_resistance)

    total_current = np.divide(voltage_series, total_resistance)
    return total_current


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

    def initialize_bypass_diode(self, ambient_celsius_temperature, irradiance_array):
        """
        Initialize the bypass diode parameters.
        """
        # Calculate the average cell temperature for the set of cells
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
        k = 1.38064852e-23  # Boltzmann constant in J/K
        q = 1.60217662e-19  # Charge of electron in Coulombs
        thermal_voltage = k * average_temperature / q
        self.bypass_diode.thermal_voltage = thermal_voltage

        # Use the first PV cell to get the saturation current and ideality factor
        first_pv_cell = self.pv_cells[0]
        # self.bypass_diode.saturation_current = first_pv_cell.reference_dark_current_density
        # self.bypass_diode.ideality_factor = first_pv_cell.gamma_ref

    def calculate_iv_curve(
        self,
        ambient_celsius_temperature: float,
        irradiance_array: np.ndarray,
        wind_speed: float,
        *,
        current_density_series: np.ndarray | None = None,
        current_series: np.ndarray | None = None,
        voltage_series: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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

        """
        self.initialize_bypass_diode(ambient_celsius_temperature, irradiance_array)

        # Calculate the curves for each cell
        cell_to_iv_series: dict[PVCell, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for pv_cell in tqdm(self.pv_cells, desc="Bypassed IV curves", leave=False):
            cell_to_iv_series[pv_cell] = pv_cell.calculate_iv_curve(
                ambient_celsius_temperature,
                irradiance_array,
                wind_speed,
                current_density_series=current_density_series,
                current_series=current_series,
                voltage_series=voltage_series,
            )
        # plt.figure()

        # Plot each cell's IV curve
        # for pv_cell in self.pv_cells:
        #     current, power, voltage = cell_to_iv_series[pv_cell]
        #     plt.plot(voltage, current, label=f'Cell {pv_cell.cell_id}')

        # plt.title('IV Curves of Individual PV Cells')
        # plt.xlabel('Voltage (V)')
        # plt.ylabel('Current (A)')
        # plt.legend()
        # plt.grid(True)
        # ylim = plt.ylim()
        # plt.show()

        # Add up the voltage for each cell
        combined_voltage_series = sum(
            cell_to_iv_series[pv_cell][2] for pv_cell in self.pv_cells
        )

        # Determine the bypass-diode curve
        bypass_diode_curve = self.bypass_diode.calculate_i_from_v(
            combined_voltage_series
        )

        # plt.plot(combined_voltage_series, cell_to_iv_series[self.pv_cells[0]][0], label='Cells', color='blue')

        # plt.plot(combined_voltage_series, bypass_diode_curve, label='Bypass Diode', color='red')

        # Calculate the total current
        total_current = bypass_diode_curve + (cell_to_iv_series[pv_cell][0])

        # Re-compute the combined power series.
        combined_power_series = total_current * combined_voltage_series

        # Plot combined IV curve (dashed)
        # plt.plot(combined_voltage_series, total_current, '--', label='Combined IV Curve')

        # Finalize plot
        # plt.xlabel('Voltage (V)')
        # plt.ylabel('Current (A)')
        # plt.legend()
        # plt.grid(True)
        # plt.title('Combined IV Curve with Bypass Diode')
        # plt.ylim(*ylim)
        # plt.show()

        import pdb

        pdb.set_trace()

        return total_current, combined_power_series, combined_voltage_series


###OLD METHOD

# from dataclasses import dataclass

# import numpy as np

# from .pv_cell import PVCell

# __all__ = ("BypassDiode", "BypassedCellString")


# @dataclass(kw_only=True)
# class BypassDiode:
#     """
#     Represents a bypass diode.

#     .. attribute:: bypass_voltage
#         The voltage at which the bypass diode will kick in and bypass the cell series.

#     .. attribute:: end_index
#         The end index for which to bypass cells.

#     .. attribute:: start_index
#         The start index for which to bypass cells.

#     """

#     bypass_voltage: float
#     end_index: int
#     start_index: int


# @dataclass(kw_only=True)
# class BypassedCellString:
#     """
#     Represents a series of cells in a string, bypassed by a single diode.

#     .. attribute:: bypass_diode
#         The bypass diode installed.

#     .. attribute:: pv_cells
#         The `list` of PV cells that are in series and bypassed by the diode.

#     """

#     bypass_diode: BypassDiode
#     pv_cells: list[PVCell]

#     @property
#     def breakdown_voltage(self) -> float:
#         """
#         Return the breakdown voltage, i.e., the bypass-diode voltage.

#         Returns:
#             - The breakdown voltage.

#         """

#         return self.bypass_diode.bypass_voltage

#     @property
#     def cell_id(self) -> float:
#         """
#         Return the cell ID for the fist cell in the string.

#         Returns:
#             - The cell ID of the first cell within the string.

#         """

#         return min([cell.cell_id for cell in self.pv_cells])

#     def __hash__(self) -> int:
#         """
#         Return a hash of the first cell ID.

#         """

#         return hash(self.cell_id)

#     def calculate_iv_curve(
#         self,
#         ambient_celsius_temperature: float,
#         irradiance_array: np.ndarray,
#         *,
#         current_density_series: np.ndarray | None = None,
#         current_series: np.ndarray | None = None,
#         voltage_series: np.ndarray | None = None,
#     ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
#         """
#         Calculate the IV curve for the bypassed string of cells.

#         Inputs:
#             - ambient_celsius_temperature:
#                 The ambient temperature, in degrees Celsius.
#             - irradiance_array:
#                 The irradiance, in W/m^2, striking all the cells in the module.
#             - current_density_series:
#                 If provided, the current-density series---the series of opints over which to
#                 calculate the current an power output from the cell.
#             - current_series:
#                 The series of current points over which to calculate the current and power
#                 output from the cell.
#             - voltage_series:
#                 The series of voltage points over which to calculate the current and power
#                 output from the cell.

#         Returns:
#             - current_series:
#                 The current values.
#             - power_series:
#                 The power values.
#             - voltage_series:
#                 The voltage series.

#         """

#         # Calculate the curves for each cell
#         cell_to_iv_series: dict[PVCell, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
#         for pv_cell in self.pv_cells:
#             cell_to_iv_series[pv_cell] = pv_cell.calculate_iv_curve(
#                 ambient_celsius_temperature,
#                 irradiance_array,
#                 current_density_series=current_density_series,
#                 current_series=current_series,
#                 voltage_series=voltage_series,
#             )

#         # Add up the voltage for each cell
#         combined_voltage_series = sum(
#             cell_to_iv_series[pv_cell][2] for pv_cell in self.pv_cells
#         )

#         # Bypass based on the diode voltage.
#         combined_voltage_series = np.array(
#             [
#                 max(entry, self.bypass_diode.bypass_voltage)
#                 for entry in combined_voltage_series
#             ]
#         )

#         # Re-compute the combined power series.
#         combined_power_series = (
#             current_series := cell_to_iv_series[self.pv_cells[0]][0]
#         ) * combined_voltage_series

#         return current_series, combined_power_series, combined_voltage_series
