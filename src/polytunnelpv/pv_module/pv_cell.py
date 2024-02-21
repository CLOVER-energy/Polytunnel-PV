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

from math import pi
from pvlib.irradiance import pvlib_get_total_irradiance


def _degree_to_radian(angle_in_degrees: float) -> float:
    """
    Converts an angle from degrees to radians.

    Inputs:
        - angle_in_degrees:
            The angle in degrees.

    Returns:
        The angle in radians.

    """

    return (pi / 180) * angle_in_degrees


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

    def __init__(
        self,
        azimuth: float,
        length: float,
        tilt: float,
        width: float,
        breakdown_voltage: float = -15,
        reference_temperature: float = 20,
    ) -> None:
        """
        Instantiate a solar cell instance.

        Inputs:
            - azimuth:
                The azimuth angle of the cell, in degrees.
            - length:
                The length of the cell, in meters.
            - tilt:
                The tilt of the cell, in degrees.
            - voltage_temperature_coefficient:
                The temperature coefficient of the open-circuit voltage.
            - width:
                The width of the cell, in meters.
            - breakdown_voltage:
                The breakdown voltage of the cell, in Volts.
            - reference_temperature:
                The reference temperature against which the cell parameters are defined,
                measured in degrees Celsius.

        """

        self.azimuth = azimuth
        self._azimuth_in_radians: float | None = None
        self.length = length
        self.tilt = tilt
        self._tilt_in_radians: float | None = None
        self.width = width
        self.breakdown_voltage = breakdown_voltage
        self.reference_temperature = reference_temperature

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
        Two cells are equal if they occupy the same position, and, hence, the same tilt.

        """

        return self.tilt == other.tilt

    def __repr__(self) -> str:
        """The default representation of the cell."""

        return self.__str__()

    def __lt__(self, other) -> bool:
        """
        The tilt can be used to determine the cell's unique identiy and sort cells.

        """

        return self.tilt < other.tilt

    def __str__(self):
        """
        A nice-looking string representation of the cell.

        Returns
            Information about the cell and light source.

        """
        return f"PVCell(azimuth={self.azimuth:.2g}, tilt={self.tilt:.2g})"

    @property
    def azimuth_in_radians(self) -> float:
        """Return the tilt in radians."""

        if self._azimuth_in_radians is None:
            self._azimuth_in_radians = _degree_to_radian(self.azimuth)

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
            self._tilt_in_radians = _degree_to_radian(self.tilt)

        return self._tilt_in_radians


def get_irradiance(
    diffuse_horizontal_irradiance: float,
    direct_normal_irradiance: float,
    pv_cell: PVCell,
    solar_azimuth: float,
    solar_zenith: float,
    global_horizontal_irradiance: float | None = None,
) -> float:

    θs, φs = (
        solar_zenith * (np.pi / 180),
        solar_azimuth * (np.pi / 180),
    )
    solar_azimuth_in_radians, solar_zenith_in_radians = map(
        _degree_to_radian, (solar_azimuth, solar_zenith)
    )
    alignment_coef = (
        np.sin(pv_cell.tilt_in_radians)
        * np.cos(pv_cell.azimuth_in_radians)
        * np.sin(θs)
        * np.cos(φs)
        + np.sin(pv_cell.tilt_in_radians)
        * np.sin(pv_cell.azimuth_in_radians)
        * np.sin(θs)
        * np.sin(φs)
        + np.cos(pv_cell.tilt_in_radians) * np.cos(θs)
    )  # using the dot product of the vector normal to the surface of the cell and the vector of the sun relative to the cell
    ghi = dhi + dni * alignment_coef
    total_irrad = pvlib_get_total_irradiance(
        pv_cell.tilt_in_radians,
        pv_cell.azimuth_in_radians,
        solar_zenith,
        solar_azimuth,
        direct_normal_irradiance,
        global_horizontal_irradiance,
        diffuse_horizontal_irradiance,
    )  # (surface_tilt, surface_azimuth, solar_zenith, solar_azimuth, dni, ghi, dhi)
    poa_global = total_irrad["poa_global"]
    return poa_global


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
