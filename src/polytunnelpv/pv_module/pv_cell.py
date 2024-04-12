#!/usr/bin/python3.10
########################################################################################
# pv_cell.py - Curved-PV-module cell Python module.                                    #
#                                                                                      #
# Author: Yaar Safra, Ben Winchester                                                   #
# Copyright: Yaar Safra, 2024                                                          #
# Date created: 21/02/2024                                                             #
# License: Open source                                                                 #
########################################################################################
"""
pv_cell.py - The pv-cell module for Polytunnel-PV.

This module provides functionality for the modelling of the individual cells within the
curved PV module.

"""

from dataclasses import dataclass
from math import cos, radians
from pvlib.irradiance import get_total_irradiance as get_total_irradiance

__all__ = ("get_irradiance", "PVCell")

# POA global key:
#   Keyword for extracting the global irradiance once computed by pvlib.
POA_GLOBAL_KEY: str = "poa_global"


@dataclass
class PVCell:
    """
    A single cell within the curved PV module.

    Attributes:
        - azimuth:
            The azimuth angle of the cell, in degrees.
        - breakdown_voltage:
            The breakdown voltage of the cell, in Volts.
        - length:
            The length of the cell, in meters.
        - reference_temperature:
            The reference temperature against which the cell parameters are defined,
            measured in degrees Celsius.
        - tilt:
            The tilt of the cell, in degrees.
        - voltage_temperature_coefficient:
            The temperature coefficient of the open-circuit voltage.
        - width:
            The width of the cell, in meters.

    """

    azimuth: float
    length: float
    tilt: float
    width: float
    breakdown_voltage: float
    reference_temperature: float = 25
    _azimuth_in_radians: float | None = None
    _cell_id: float | int | None = None
    _tilt_in_radians: float | None = None

    # c_params = {
    #     "I_L_ref": 8.24,
    #     "I_o_ref": 2.36e-9,
    #     "a_ref": 1.3 * Vth,
    #     "R_sh_ref": 1000,
    #     "R_s": 0.00181,
    #     "alpha_sc": 0.0042,
    #     "breakdown_factor": 2e-3,
    #     "breakdown_exp": 3,
    #     "breakdown_voltage": self._breakdown_voltage,
    # }
    # self._c_params = c_params

    def __eq__(self, other) -> bool:
        """
        Two cells are equal if they occupy the same space, i.e., are oriented the same.

        """

        return (self.tilt == float(other.tilt)) and (
            self.azimuth == float(other.azimuth)
        )

    def __repr__(self) -> str:
        """The default representation of the cell."""

        return self.__str__()

    @property
    def cell_id(self) -> float | int:
        """A unique ID for the cell, based on the orientation if not provided."""

        if self._cell_id is None:
            self._cell_id = self.tilt

        return self._cell_id

    def __hash__(self) -> int:
        """Return a unique number based on the cell ID."""

        return hash(self.cell_id)

    def __lt__(self, other) -> bool:
        """
        The tilt, combined with the azimuth, can be used to determine the cell's unique
        identiy and sort cells.

        """

        return self.cell_id < float(other.cell_id)

    def __str__(self) -> str:
        """
        A nice-looking string representation of the cell.

        Returns
            Information about the cell and light source.

        """
        return f"PVCell(azimuth={self.azimuth:.2f}, tilt={self.tilt:.2f})"

    @property
    def azimuth_in_radians(self) -> float:
        """Return the tilt in radians."""

        if self._azimuth_in_radians is None:
            self._azimuth_in_radians = radians(self.azimuth)

        return self._azimuth_in_radians

    def bypass(self, bypass_diode_voltage: float | None = None):
        if bypass_diode_voltage is None:
            raise Exception(
                "Must call `bypass` function with a new bypass-diode voltage."
            )

        self.set_breakdown_voltage(bypass_diode_voltage)

    def set_breakdown_voltage(self, breakdown_voltage_to_set) -> None:
        self.breakdown_voltage = breakdown_voltage_to_set

    @property
    def tilt_in_radians(self) -> float:
        """Return the tilt in radians."""

        if self._tilt_in_radians is None:
            self._tilt_in_radians = radians(self.tilt)

        return self._tilt_in_radians


def get_irradiance(
    pv_cell: PVCell,
    diffuse_horizontal_irradiance: float,
    global_horizontal_irradiance: float,
    solar_azimuth: float,
    solar_zenith: float,
    direct_normal_irradiance: float | None = None,
) -> float | None:
    """
    Compute the irradiance falling on an individual solar cell.

    Inputs:
        - pv_cell:
            The cell to calculate these values for.
        - diffuse_horizontal_irradiance:
            The diffuse horizontal irradiance in W/m^2.
        - global_horizontal_irradiance:
            The global horizontal irradiance.
        - solar_azimuth:
            The current azimuth angle of the sun.
        - solar_zenith:
            The current zenith angle of the sun, _i.e._, it's declanation.
        - direct_normal_irradiance:
            If provided, the direct normal irradiance in W/m^2. If not, this is
            calculated using the solar zenith angle.

    """

    # If it's nighttime, _i.e._, the sun is below the horizon or the global irradiance
    # is zero, then simply return zero.
    if solar_zenith >= 90 or global_horizontal_irradiance <= 0:
        return 0

    # Determine the DNI from the GHI and DHI.
    if direct_normal_irradiance is None:
        direct_normal_irradiance = (
            global_horizontal_irradiance - diffuse_horizontal_irradiance
        ) / cos(radians(solar_zenith))

    # Call to PVlib to calculate the total irradiance incident on the surface.
    total_irradiance = get_total_irradiance(
        pv_cell.tilt,
        pv_cell.azimuth,
        solar_zenith,
        solar_azimuth,
        direct_normal_irradiance,
        global_horizontal_irradiance,
        diffuse_horizontal_irradiance,
    )

    # Extract and return the global irradiance striking the surface.
    return total_irradiance.get(POA_GLOBAL_KEY, None)  # type: ignore [no-any-return]


def get_iv_curve(self, show_axis=True):
    curves = []
    labels = []
    for i in range(len(SC._cell_list)):
        curve = simulate_full_curve(
            SC._cell_list[i].get_c_params(),
            SC.get_irradiance(SC._cell_list[i]),
            SC._cell_list[i]._temp,
            break_volt=SC._cell_list[i].get_breakdown_voltage(),
        )
        curves.append(curve)
        labels.append("angle =" + str(round(SC.get_cell_tilt(SC._cell_list[i]), 2)))
    axis = plot_curves(curves, labels)
    if show_axis:
        return axis
    else:
        return curves
