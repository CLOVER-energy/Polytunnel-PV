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

import enum
import warnings

from dataclasses import dataclass
from math import radians

import numpy as np

import pvlib


__all__ = ("CellType", "get_irradiance", "PVCell", "relabel_cell_electrical_parameters")

# A_REF:
#   Keyword for the a-ref parameter.
A_REF: str = "a_ref"

# DIODE_IDEALITY_FACTOR:
#   Keyword for the diode ideality factor.
DIODE_IDEALITY_FACTOR: str = "_reference_diode_ideality_factor"

# NUM_CELLS_IN_PARENT_MODULE:
#   Keyword for the number of cells in the parent module.
NUM_CELLS_IN_PARENT_MODULE: str = "_num_cells_in_parent_module"

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

# STEFAN_BOLTZMAN_CONSTANT:
#   The Stefan-Boltzman constant in SI units.
STEFAN_BOLTZMAN_CONSTANT: float = 5.670374419 * (10**-8)

# TEMPERATURE_PRECISION:
#   The precision required when solving the matrix equation for the system temperatures.
TEMPERATURE_PRECISION: float = 0.1

# ZERO_CELSIUS_OFFSET:
#   The offset to use for converting to Kelvin.
ZERO_CELSIUS_OFFSET: float = 273.15


def _conductive_air_heat_transfer_coefficient(wind_speed: float) -> float:
    """
    Calculates the conductive heat-transfer coefficient between PV panels and the air.

    The heat transfer coefficient is given by:
        h_air = 3.97(13) sWm^3K v_wind
              + 6.90(5) W/m^2K

    """

    return float(3.97 * wind_speed + 6.90)


def _radiation_to_sky_coefficient(
    collector_emissivity: float, collector_temperature: float, sky_temperature: float
) -> float:
    """
    Returns the heat-transfer coefficient between a collector and the sky in W/m^2.

    This function has been taken, and reproduced with permission, from
    Winchester, B., Nelson, J., and Markides, C.N., "HEATDesalination" [Software]
    Available from: https://github.com/BenWinchester/HEATDesalination

    :param: **collector_emissivity:**
        The emissivity of the collector.

    :param: **collector_temperature:**
        The temperature of the collector, measured in Kelvin.

    :param: **sky_temperature:**
        The sky temperature, measured in Kelvin.

    :returns:
        The heat-transfer coefficient between the collector and the sky, measured in
        Watts per meter squared.

    """
    return (
        STEFAN_BOLTZMAN_CONSTANT  # [W/m^2*K^4]
        * collector_emissivity
        * (collector_temperature**2 + sky_temperature**2)  # [K^2]
        * (collector_temperature + sky_temperature)  # [K]
    )


def _sky_temperature(ambient_temperature: float) -> float:
    """
    Determines the radiative temperature of the sky.

    This function has been taken, and reproduced with permission, from
    Winchester, B., Nelson, J., and Markides, C.N., "HEATDesalination" [Software]
    Available from: https://github.com/BenWinchester/HEATDesalination

    The "sky," as a black body, has a radiative temperature different to that of the
    surrounding air, or the ambient temperature. This function converts between them and
    outputs the sky's radiative temperature.

    :returns: The radiative temperature of the "sky" in Kelvin.

    """

    return float(0.0552 * (ambient_temperature**1.5))


class CellType(enum.Enum):
    """
    Used to denote whether the cells are mono- or bi-facial in nature.

    - MONO_FACIAL: 0
        Used to denote mono-facial cells.

    - BIFACIAL: 1
        Used to denote bi-facial cells.

    """

    MONO_FACIAL: int = 0
    BIFACIAL: int = 1


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

    # Private attributes:
    #
    # .. attribute:: _azimuth_in_radians:
    #   The azimuth angle in radians.
    #
    # .. attribute:: _reference_diode_ideality_factor:
    #   The ideality factor for the diode.
    #
    # .. attribute:: _cell_id:
    #   The ID number for the cell.
    #
    # .. attribute:: num_cells_in_parent_module:
    #   The number of cells in the parent module from which the electrical parameters
    #   for the cell are obtained.
    #   NOTE: This is used purely in electrical calculations and should not be used
    #   otherwise.
    #
    # .. attribute:: _tilt_in_radians:
    #   The tilt angle in radians.
    #
    # .. attribute:: __open_circuit_voltage:
    #   The open-circuit voltage, in Volts.
    #
    # .. attribute:: __mpp_thermal_coefficient:
    #   The thermal coefficient of the maximum power point.
    #
    # .. attribute:: __reference_efficiency:
    #   The power-conversion efficiency of the cell at reference conditions.
    #

    azimuth: float
    cell_type: CellType
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
    absorptivity: float = 0.9
    emissivity: float = 0.9
    _azimuth_in_radians: float | None = None
    _cell_id: float | int | None = None
    _reference_diode_ideality_factor: float | int | None = None
    _num_cells_in_parent_module: int | None = None
    _tilt_in_radians: float | None = None
    __open_circuit_voltage: float | None = None
    __mpp_thermal_coefficient: float | None = None
    __reference_efficiency: float | None = None
    __short_circuit_current: float | None = None

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
    def num_cells_in_parent_module(self) -> float:
        """Return the number of cells in the parent module."""

        if self._num_cells_in_parent_module is None:
            return 1

        return self._num_cells_in_parent_module

    @property
    def alpha_sc(self) -> float:
        """The short-circuit curent density tempearture coefficient."""

        return self.short_circuit_current_density_temperature_coefficient

    @property
    def area(self) -> float:
        """The area of the cell in meters squared."""

        return self.width * self.length

    @property
    def gamma_ref(self) -> float:
        """The reference diode ideality factor."""

        if self._reference_diode_ideality_factor is None:
            return 1

        return self._reference_diode_ideality_factor

    @property
    def mpp_thermal_coefficient(self) -> float:
        """
        Calculate the thermal coefficient of the maximum power point.

        Returns:
            The maximum power point of the thermal coefficient.

        """

        # Calculate the MPP thermal coefficient if it hasn't been calculated yet.
        if self.__mpp_thermal_coefficient is None:
            self.__mpp_thermal_coefficient = pvlib.ivtools.sdm.pvsyst_temperature_coeff(
                self.alpha_sc,
                self.gamma_ref,
                0,
                self.reference_photocurrent_density,
                self.reference_dark_current_density,
                self.reference_shunt_resistance,
                self.reference_shunt_resistance,
                self.reference_series_resistance,
                self.num_cells_in_parent_module,
                EgRef=self.eg_ref,
                irrad_ref=self.reference_irradiance,
                temp_ref=self.reference_temperature,
            )

        return self.__mpp_thermal_coefficient

    @property
    def open_circuit_voltage(self) -> float:
        """
        Calculate (if necessary) and return the open-circuit voltage of the cell.

        Returns:
            The open-circuit voltage of the cell.

        """

        if self.__open_circuit_voltage is None:
            _calculated_pv_cell_params = pvlib.pvsystem.calcparams_desoto(
                self.reference_irradiance,
                self.reference_temperature,
                self.alpha_sc,
                self.a_ref,
                self.j_l_ref,
                self.j_o_ref,
                self.r_sh_ref,
                self.r_s,
                self.eg_ref,
                self.d_eg_dt_ref,
            )

            reference_open_circuit_voltage = pvlib.pvsystem.singlediode(
                *_calculated_pv_cell_params
            )["v_oc"]

            self.__open_circuit_voltage = reference_open_circuit_voltage

        return self.__open_circuit_voltage

    @property
    def reference_efficiency(self) -> float:
        """
        Calculate an approximate reference efficiency value.

        Returns:
            An approximate value of the reference efficiency of the cell.

        """

        if self.__reference_efficiency is None:
            # Compute the efficiency of the cell at reference conditions.
            _calculated_pv_cell_params = pvlib.pvsystem.calcparams_desoto(
                self.reference_irradiance,
                self.reference_temperature,
                self.alpha_sc,
                self.a_ref,
                self.j_l_ref,
                self.j_o_ref,
                self.r_sh_ref,
                self.r_s,
                self.eg_ref,
                self.d_eg_dt_ref,
            )

            _, _, maximum_power_point_power = pvlib.singlediode.bishop88_mpp(
                *_calculated_pv_cell_params
            )

            self.__reference_efficiency = maximum_power_point_power / (
                self.reference_irradiance * self.area * self.num_cells_in_parent_module
            )

        return self.__reference_efficiency

    @property
    def short_circuit_current(self) -> float:
        """
        Calculate (if necessary) and return the short-circuit current of the cell.

        Returns:
            The short-circuit current of the cell in Amps.

        """

        if self.__short_circuit_current is None:
            _calculated_pv_cell_params = pvlib.pvsystem.calcparams_desoto(
                self.reference_irradiance,
                self.reference_temperature,
                self.alpha_sc,
                self.a_ref,
                self.j_l_ref,
                self.j_o_ref,
                self.r_sh_ref,
                self.r_s,
                self.eg_ref,
                self.d_eg_dt_ref,
            )

            reference_short_circuit_current = pvlib.pvsystem.singlediode(
                *_calculated_pv_cell_params
            )["i_sc"]

            self.__short_circuit_current = reference_short_circuit_current

        return self.__short_circuit_current

    @property
    def short_circuit_current_density(self) -> float:
        """
        Return the short-circuit current density of the cell in Amps per meter squared.

        """

        return self.short_circuit_current / self.area

    def average_cell_temperature(
        self,
        ambient_temperature: float,
        solar_irradiance: float,
        wind_speed: float,
    ) -> float:
        """
        Calculate the temperature of the PV module in Kelvin.

        This function has been taken, and reproduced with permission, from
        Winchester, B., Nelson, J., and Markides, C.N., "HEATDesalination" [Software]
        Available from: https://github.com/BenWinchester/HEATDesalination

        Uses a heat-balance calculation iteratively to solve for the average temperature of
        the PV module:
            0W = a_pv G (1 - eta_el)
                 + h_air (T_air - T_pv)
                 + e' sigma (T_sky - Tpv),
            where:
                a_pv    is the absorptivity of the collector,
                G       the solar irradiance in Watts per meter squared,
                eta_el  the electrical efficiency of the collector,
                h_air   the conductive heat-transfer coefficient between the collector
                        and the surrounding air,
                e'      the effective emissivity made linear by combining temperature
                        terms of higher orders,
            and sigma   is the Stefan-Boltzman coefficient.

        This can be rearranged to the form employed here:
            T_pv = (
                a_pv G (1 - eta_ref (1 + beta T_ref))
                + h_air T_amb
                + e' sigma T_sky
            ) / (
                e' sigma + h_air - a_pv beta eta_ref G
            )
        which is linear in temperature and can be solved iteratively as the value of e'
        contains higher-order terms in temperature.

        :param: **ambient_temperature:**
            The ambient temperature, measured in degrees Kelvin.

        :param: **solar_irradiance:**
            The solar irradiance incident on the collector, measured in Watts per meter
            squared.

        :param: **wind_speed:**
            The wind speed in meters per second.

        :returns:
            The average temperature of the PV module installed in degrees Kelvin.

        """

        # Calculate variable values which remain constant throughout the iteration
        sky_temperature: float = _sky_temperature(ambient_temperature)

        # Setup inputs for the iterative loop
        best_guess_average_temperature: float = ambient_temperature
        solution_found: bool = False

        # Loop through until a solution is found
        while not solution_found:
            # Compute the necessary coefficients
            radiation_to_sky_coefficient = _radiation_to_sky_coefficient(
                self.emissivity, best_guess_average_temperature, sky_temperature
            )

            # Calculate the average temperature of the collector
            average_temperature: float = (
                (self.absorptivity * solar_irradiance)
                * (
                    1
                    - self.reference_efficiency
                    * (1 + self.mpp_thermal_coefficient * self.reference_temperature)
                )
                + _conductive_air_heat_transfer_coefficient(wind_speed)
                * ambient_temperature
                + radiation_to_sky_coefficient * sky_temperature
            ) / (
                radiation_to_sky_coefficient
                + _conductive_air_heat_transfer_coefficient(wind_speed)
                - self.absorptivity
                * self.mpp_thermal_coefficient
                * self.reference_efficiency
                * solar_irradiance
            )

            # Break if this average temperature is within the required precision.
            if (
                abs(best_guess_average_temperature - average_temperature)
                < TEMPERATURE_PRECISION
            ):
                solution_found = True

            best_guess_average_temperature = average_temperature

        return average_temperature

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

    def rescale_voltage(
        self, voltage_to_rescale: float | list[float] | np.ndarray
    ) -> float | np.ndarray:
        """Rescale voltage values based on the cell being in series."""

        def _rescale_voltage(voltage: float) -> float:
            return voltage / self.num_cells_in_parent_module

        if isinstance(voltage_to_rescale, float):
            return _rescale_voltage(voltage_to_rescale)

        return np.asarray([_rescale_voltage(voltage) for voltage in voltage_to_rescale])

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
        Calculate the IV curve for the cell.

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
            If the cell has broken down, this is returned.

        """

        cell_temperature = (
            self.average_cell_temperature(
                ambient_celsius_temperature + ZERO_CELSIUS_OFFSET,
                (solar_irradiance := irradiance_array.iloc[self.cell_id]),
                wind_speed,
            )
            - ZERO_CELSIUS_OFFSET
        )

        current_series, voltage_series, power_series = calculate_cell_iv_curve(
            cell_temperature,
            solar_irradiance,
            self,
            current_density_series=current_density_series,
            current_series=current_series,
            voltage_series=voltage_series,
        )

        return (
            current_series,
            voltage_series,
            power_series,
            list(voltage_series >= self.breakdown_voltage),
        )


def calculate_cell_iv_curve(
    cell_temperature: float,
    irradiance: float,
    pv_cell: PVCell,
    *,
    current_density_series: np.ndarray | None = None,
    current_series: np.ndarray | None = None,
    voltage_series: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the IV curve for a cell.

    NOTE: Depending on the units that the user uses, this function can return a result
    with units of current OR current density. It _should_ be implemented with units of
    current density, so a user should be passing in reference J values.

    NOTE: The current OR the voltage series should be supplied, not both, and this
    function will calculate whichever is not supplied.

    :param: **cell_temperature:**
        The cell temperature, in degrees Celsius.

    :param: **irradiance:**
        The irradiance, in W/m^2, striking the cell.

    :param: **pv_cell:**
        The PV cell being considered.

    :param: **current_density_series:**
        If provided, the current-density series---the series of opints over which to
        calculate the current an power output from the cell.

    :param: **current_series:**
        The series of current points over which to calculate the current and power
        output from the cell.

    :param: **voltage_series:**
        The series of voltage points over which to calculate the current and power
        output from the cell.

    :returns: current_series
        The current values.

    :returns: power_series:
        The power values.

    :returns: voltage_series:
        The voltage series.

    """

    if (
        current_series is not None or current_density_series is not None
    ) and voltage_series is not None:
        raise Exception(
            "Must supply either the current or voltage ranges to calculate over, not "
            "both, as this is ambiguous."
        )

    if (
        current_series is None and current_density_series is None
    ) and voltage_series is None:
        raise Exception(
            "Must supply either the current or voltage ranges to calculate over: "
            "neither has been supplied."
        )

    _calculated_pv_cell_params = pvlib.pvsystem.calcparams_desoto(
        irradiance,
        cell_temperature,
        pv_cell.alpha_sc,
        pv_cell.a_ref,
        pv_cell.j_l_ref,
        pv_cell.j_o_ref,
        pv_cell.r_sh_ref,
        pv_cell.r_s,
        pv_cell.eg_ref,
        pv_cell.d_eg_dt_ref,
    )

    def bishop88_current_wrapper_function(*args, **kwargs):
        """
        Function to wrap around the bishop 88 method.

        The `pvlib.singlediode.bishop88_i_from_v` method throws `RuntimeError`s when
        the current would result in a value that is out of bounds.
        Rather than throw this error to the user, it is caught here, and infinities
        are returned instead.

        Return:
            - The result of the call to `pvlib.singlediode.bishop88_i_from_v`,
                provided that it's an allowed results;
            - `np.inf` or `-np.inf` otherwise.

        """

        with warnings.catch_warnings():
            warnings.filterwarnings("error")
            try:
                current = pvlib.singlediode.bishop88_i_from_v(*args, **kwargs)
            except (RuntimeError, RuntimeWarning):
                return np.inf
            else:
                if current > 0:
                    return current
                return -np.inf

    def bishop_voltage_wrapper_function(*args, **kwargs):
        """
        Function to wrap around the bishop 88 method.

        The `pvlib.singlediode.bishop88_v_from_i` method throws `RuntimeError`s when
        the current would result in a value that is out of bounds.
        Rather than throw this error to the user, it is caught here, and infinities
        are returned instead.

        Return:
            - The result of the call to `pvlib.singlediode.bishop88_v_from_i`,
                provided that it's an allowed results;
            - `np.inf` or `-np.inf` otherwise.

        """

        with warnings.catch_warnings():
            warnings.filterwarnings("error")
            try:
                voltage = pvlib.singlediode.bishop88_v_from_i(*args, **kwargs)
            except (RuntimeError, RuntimeWarning):
                return 0
            else:
                if voltage > (
                    module_breakdown_voltage := (
                        pv_cell.num_cells_in_parent_module * pv_cell.breakdown_voltage
                    )
                ):
                    return voltage
                return module_breakdown_voltage

    if voltage_series is None and (
        current_density_series is not None or current_series is not None
    ):
        # Use the current density if provided, otherwise the current
        voltage_series = pv_cell.rescale_voltage(  # type: ignore [assignment]
            [
                bishop_voltage_wrapper_function(
                    current,
                    *_calculated_pv_cell_params,
                    breakdown_voltage=pv_cell.breakdown_voltage,
                    breakdown_factor=2e-3,
                    breakdown_exp=3,
                )
                for current in (
                    current_series  # type: ignore [union-attr]
                    if current_series is not None
                    else current_density_series
                )
            ]
        )
    elif voltage_series is not None:
        current_series = np.asarray(
            [
                bishop88_current_wrapper_function(
                    voltage,
                    *_calculated_pv_cell_params,
                    breakdown_voltage=pv_cell.breakdown_voltage,
                    breakdown_factor=2e-3,
                    breakdown_exp=3,
                )
                for voltage in voltage_series
            ]
        )
    else:
        raise Exception(
            "Must specify either a current-density, current, or voltage series."
        )

    if current_density_series is not None:
        current_series = current_density_series * pv_cell.area

    if current_series is None:
        raise Exception(
            "An internal error occurred when computing the current series. Debug."
        )

    if voltage_series is None:
        raise Exception(
            "An internal error occurred when computing the voltage series. Debug."
        )

    power_series = current_series * voltage_series

    return current_series, power_series, voltage_series


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

    :param: **pv_cell:
        The cell to calculate these values for.

    :param: **diffuse_horizontal_irradiance:
        The diffuse horizontal irradiance in W/m^2.

    :param: **global_horizontal_irradiance:
        The global horizontal irradiance.

    :param: **solar_azimuth:
        The current azimuth angle of the sun.

    :param: **solar_zenith:
        The current zenith angle of the sun, _i.e._, it's declanation.

    :param: **direct_normal_irradiance:
        If provided, the direct normal irradiance in W/m^2. If not, this is calculated
        using the solar zenith angle.

    :returns: The irradiance striking the cell.

    """

    # If it's nighttime, _i.e._, the sun is below the horizon or the global irradiance
    # is zero, then simply return zero.
    if solar_zenith >= 90 or global_horizontal_irradiance <= 0:
        return 0

    # Determine the pvlib.irradiance.dni from the GHI and DHI.
    if direct_normal_irradiance is None:
        raise Exception("Internal error encountered in irradiance calculation.")
        # direct_normal_irradiance = pvlib.irradiance.dni(
        #     global_horizontal_irradiance, diffuse_horizontal_irradiance, solar_zenith
        # )

    # Call to PVlib to calculate the total irradiance incident on the surface.
    total_irradiance = pvlib.irradiance.get_total_irradiance(
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
    cell_electrical_params: dict[str, float],
) -> dict[str, float]:
    """
    Re-map entries within the mapping to more-sensible and user-readable variable names.

    :param: The electrical parameters governing a pv cell.

    :returns:
        The electrical parameters governing a pv cell, with variable names processed.

    """

    return {
        A_REF: cell_electrical_params["a_ref"],
        DIODE_IDEALITY_FACTOR: cell_electrical_params["gamma_r"],
        NUM_CELLS_IN_PARENT_MODULE: cell_electrical_params["N_s"],
        REFERENCE_DARK_CURRENT_DENSITY: cell_electrical_params["I_o_ref"],
        # / cell_electrical_params["A_c"],
        REFERENCE_PHOTOCURRENT_DENSITY: cell_electrical_params["I_L_ref"],
        # / cell_electrical_params["A_c"],
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
