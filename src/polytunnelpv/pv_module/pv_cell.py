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
from math import radians
from pvlib.irradiance import dni, get_total_irradiance

__all__ = ("get_irradiance", "PVCell", "relabel_cell_electrical_parameters")

# A_REF:
#   Keyword for the a-ref parameter.
A_REF: str = "a_ref"

# POA global key:
#   Keyword for extracting the global irradiance once computed by pvlib.
POA_GLOBAL_KEY: str = "poa_global"

# REFERENCE_DARK_CURRENT_DENSITY:
#   Keyword for the reference dark-current-density parameter.
REFERENCE_DARK_CURRENT_DENSITY: str = "reference_dark_current_density"

# REFERENCE_PHOTOCURRENT_DENSITY:
#   Keyword for the reference photocurrent density.
REFERENCE_PHOTOCURRENT_DENSITY: str = "reference_photocurrent_density"

# REFERENCE_SERIES_RESISTANCE:
#   Keyword for the reference series resistance.
REFERENCE_SERIES_RESISTANCE: str = "reference_series_resistance"

# REFERENCE_SHUNT_RESISTANCE:
#   Keyword for the reference series resistance.
REFERENCE_SHUNT_RESISTANCE: str = "reference_shunt_resistance"

# SHORT_CIRCUIT_CURRENT_DENSITY_TEMPERATURE_COEFFICIENT:
#   Keyword for the short-circuit current-density temperature coefficient.
SHORT_CIRCUIT_CURRENT_DENSITY_TEMPERATURE_COEFFICIENT: str = (
    "short_circuit_current_density_temperature_coefficient"
)


@dataclass(kw_only=True)
class PVCell:
    """
    A single cell within the curved PV module.

    Parameters that govern the IV (current-voltage) curve of the cell, as required to
    calculate these curves, are included as attributes of the cell as they will depend
    on the cell material used.

    Attributes:
        - azimuth:
            The azimuth angle of the cell, in degrees.
        - length:
            The length of the cell, in meters.
        - tilt:
            The tilt of the cell, in degrees.
        - width:
            The width of the cell, in meters.
        - breakdown_voltage:
            The breakdown voltage of the cell, in Volts.

    IV-curve-related attributes:
        NOTE: These are taken, for the most part, from the descriptions provided in the
        pvlib library.
        - a_ref:
            The product between the diode ideality factor and the cell's thermal
            voltage.
            NOTE: This value should be cell-specific, not module-wide.
        - reference_dark_current_density:
            The reference dark (or reverse-saturation) curren, measured in Amps per
            meter squared.
        - reference_photocurrent_density:
            The reference photocurrent density, measured in Amps per meter squared.
        - reference_series_resistance:
            The series resistance at reference conditions, in ohms.
        - reference_shunt_resistance:
            The shunt resistance at reference conditions, in ohms.
        - short_circuit_current_density_temperature_coefficient:
            The temperature coefficient of the short-circuit current, in Amps per
            degree.
        - reference_bandgap_energy:
            The energy of the bandgap at reference conditions, measured in eV.
            NOTE: 1.121 eV is the default value for crystalline Silicon.
        - reference_bandgap_energy_temperature_coefficient:
            The temperature dependence of the energy bandgap at reference conditions.
            NOTE: -0.0002677 is the default value for cyrstalline Silicon.
        - reference_irradiance:
            The irradiance at reference conditions, measured in Watts per meter squared.
        - reference_temperature:
            The cell temeprature at reference conditions, measured in Watts per meter
            squared.

    """

    azimuth: float
    length: float
    tilt: float
    width: float
    breakdown_voltage: float
    a_ref: float
    reference_dark_current_density: float
    reference_photocurrent_density: float
    reference_series_resistance: float
    reference_shunt_resistance: float
    short_circuit_current_density_temperature_coefficient: float
    reference_bandgap_energy: float = 1.121
    reference_bandgap_energy_temperature_coefficient: float = -0.0002677
    reference_irradiance: float = 1000
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
    def alpha_sc(self) -> float:
        """The short-circuit curent density tempearture coefficient."""

        return self.short_circuit_current_density_temperature_coefficient

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

    @property
    def d_eg_dt_ref(self) -> float:
        """The reference temperature dependence of the bandgap energy."""

        return self.reference_bandgap_energy_temperature_coefficient

    @property
    def eg_ref(self) -> float:
        """The reference temperature dependence of the bandgap energy."""

        return self.reference_bandgap_energy_temperature_coefficient

    @property
    def j_l_ref(self) -> float:
        """The reference photocurrent density."""

        return self.reference_photocurrent_density

    @property
    def j_o_ref(self) -> float:
        """The reference dark-current density."""

        return self.reference_dark_current_density

    @property
    def r_s(self) -> float:
        """The reference series resistance."""

        return self.reference_series_resistance

    @property
    def r_sh_ref(self) -> float:
        """The reference shunt resistance."""

        return self.reference_shunt_resistance

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
        direct_normal_irradiance = dni(
            global_horizontal_irradiance, diffuse_horizontal_irradiance, solar_zenith
        )

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


def relabel_cell_electrical_parameters(
    cell_electrical_params: dict[str, float]
) -> dict[str, float]:
    """
    Re-map entries within the mapping to more-sensible and user-readable variable names.

    Inputs:
        The electrical parameters governing a pv cell.

    Returns:
        The electrical parameters governing a pv cell, with variable names processed.

    """

    return {
        A_REF: cell_electrical_params["a_ref"],  # / cell_electrical_params["N_s"],
        REFERENCE_DARK_CURRENT_DENSITY: cell_electrical_params["I_o_ref"]
        / cell_electrical_params["A_c"],
        REFERENCE_PHOTOCURRENT_DENSITY: cell_electrical_params["I_L_ref"]
        / cell_electrical_params["A_c"],
        REFERENCE_SERIES_RESISTANCE: cell_electrical_params["R_s"],
        REFERENCE_SHUNT_RESISTANCE: cell_electrical_params["R_sh_ref"],
        SHORT_CIRCUIT_CURRENT_DENSITY_TEMPERATURE_COEFFICIENT: cell_electrical_params[
            "alpha_sc"
        ],
    }


# def get_iv_curve(self, show_axis=True):
#     curves = []
#     labels = []
#     for i in range(len(SC._cell_list)):
#         curve = (
#             SC._cell_list[i].get_c_params(),
#             SC.get_irradiance(SC._cell_list[i]),
#             SC._cell_list[i]._temp,
#             break_volt=SC._cell_list[i].get_breakdown_voltage(),
#         )
#         curves.append(curve)
#         labels.append("angle =" + str(round(SC.get_cell_tilt(SC._cell_list[i]), 2)))
#     axis = plot_curves(curves, labels)
#     if show_axis:
#         return axis
#     else:
#         return curves
