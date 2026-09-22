########################################################################################
# msi_pv.py --- Module for computing the performance of m-Si PV modules.               #
#                                                                                      #
# Author(s): Benedict Winchester                                                       #
# Date created: Summer 2026                                                            #
#                                                                                      #
########################################################################################

"""
msi_pv.py - The m-Si module for the Polytunnel-PV model.

The model computers the output power from the curved modules on a polytunnel. This
module is responsible for computing the output power from conventional flat m-Si modules
and is taken with permission from HEATDesalation.

"""

import abc
import dataclasses
import enum

from logging import Logger
from typing import Any, Tuple, Type

__all__ = (
    "PerformanceCurve",
    "PVPanel",
)

# AREA:
#   Keyword for the area of the panel.
AREA: str = "area"

# COST:
#   Keyword for the cost of a component.
COST: str = "cost"

# ELECTRIC_PERFORMANCE_CURVE:
#   Keyword for the electric performance curve.
ELECTRIC_PERFORMANCE_CURVE: str = "electric_performance_curve"

# EMISSIONS:
#   Keyword for the emissions of a component.
EMISSIONS: str = "emissions"

# EMISSIONS_RANGE:
#   Keyword for the uncertainty range associated with the emissions of a component.
EMISSIONS_RANGE: str = "emissions_range"

# FIRST_ORDER:
#   Keyword for parsing first-order coefficient.
FIRST_ORDER: str = "first_order"

# LAND_USE:
#   Keyword for the land use of each panel.
LAND_USE: str = "land_use"

# MAX_MASS_FLOW_RATE:
#   Keyword for the maximum mass flow rate of HTF through the panel.
MAX_MASS_FLOW_RATE: str = "max_mass_flow_rate"

# MIN_MASS_FLOW_RATE:
#   Keyword for the minimum mass flow rate of HTF through the panel.
MIN_MASS_FLOW_RATE: str = "min_mass_flow_rate"

# NAME:
#   Keyword for parsing the name of the object.
NAME: str = "name"

# NOMINAL_MASS_FLOW_RATE:
#   Keyword for the nominal mass flow rate of HTF through the panel.
NOMINAL_MASS_FLOW_RATE: str = "nominal_mass_flow_rate"

# NOMINAL_POWER:
#   Keyword for the nominal power of the panel.
NOMINAL_POWER: str = "nominal_power"

# PV_MODULE_CHARACTERISTICS:
#   Keyword for module characteristics.
PV_MODULE_CHARACTERISTICS: str = "pv_module_characteristics"

# PV_UNIT:
#   Keyword for the unit of power outputted per PV panel.
PV_UNIT: str = "pv_unit"

# REFERENCE_EFFICIENCY:
#   Keyword for the reference efficiency of a PV panel.
REFERENCE_EFFICIENCY: str = "reference_efficiency"

# REFERENCE_SOLAR_IRRADIANCE:
#   The reference solar irradiance at which collectors are rated.
REFERENCE_SOLAR_IRRADIANCE: float = 1000

# REFERENCE_TEMPERATURE:
#   Keyword for the reference temperature of a PV panel.
REFERENCE_TEMPERATURE: str = "reference_temperature"

# SECOND_ORDER:
#   Keyword for parsing second-order coefficient.
SECOND_ORDER: str = "second_order"

# STAGNATION_TEMPERATURE:
#   Keyword for the stagnation temperature of the collector.
STAGNATION_TEMPERATURE: str = "stagnation_temperature"

# STEFAN_BOLTZMAN_CONSTANT:
#   The Stefan-Boltzman constant in SI units.
STEFAN_BOLTZMAN_CONSTANT: float = 5.670374419 * (10**-8)

# TEMPERATURE_PRECISION:
#   The precision required when solving the matrix equation for the system temperatures.
TEMPERATURE_PRECISION: float = 0.1

# THERMAL_COEFFICIENT:
#   Keyword for the temperature coefficient for the performance of a PV panel.
THERMAL_COEFFICIENT: str = "thermal_coefficient"

# THERMALPERFORMANCE_CURVE:
#   Keyword for the thermal performance curve.
THERMAL_PERFORMANCE_CURVE: str = "thermal_performance_curve"

# ZERO_CELCIUS_OFFSET:
#   Keyword for the offset of Kelvin to Celcius.
ZERO_CELCIUS_OFFSET: float = 273.15

# ZEROTH_ORDER:
#   Keyword for parsing zeroth-order coefficient.
ZEROTH_ORDER: str = "zeroth_order"


class CostableComponent:
    """
    Represents a costable component.

    NOTE: As this is a base class, public methods are not defined here.

    .. attribute:: cost
        The cost of the component, per unit component.

    .. attribute:: name
        The name of the component.

    """

    def __init__(self, cost: float, name: str) -> None:
        """
        Instantiate the costable component.

        Inputs:
            - cost:
                The cost of the component, per unit component.
            - name:
                The name of the component.

        """

        self.cost = cost
        self.name = name


class EmissableComponent(CostableComponent):
    """
    Represents an emissable component.

    .. attribute:: emissions
        The carbon emissions associated with this component, per unit component, arising
        from the manufacturing (and transport), i.e., embedded/embodied in the
        component measured in kgCO2-eq.

    .. attribute:: emissions_range
        The range of uncertainty in the carbon emissions associated witht the component,
        measured in kgCO2-eq.

    """

    def __init__(
        self, cost: float, emissions: float, name: str, *, emissions_range: float = 0
    ) -> None:
        """
        Instantiate the costable component.

        Inputs:
            - cost:
                The cost of the component, per unit component.
            - emissions:
                The carbon emissions associated with the component.
            - emissions_range:
                The range of uncertainty in the carbon emissions associated with the
                component.
            - name:
                The name of the component.

        """

        super().__init__(cost, name)

        self.emissions = emissions
        self.emissions_range = emissions_range


class ProgrammerJudgementFault(Exception):
    """
    Raised when an error occurs in the code.

    .. attribute:: code_location
        The location within the code where the error occurred.

    """

    def __init__(self, code_location: str, msg: str) -> None:
        """
        Instantiate a :class:`ProgrammerJudgementFault`.

        Inputs:
            - code_location:
                The location within the code where the error occurred.
            - msg:
                The message to append.

        """

        self.code_location = code_location
        super().__init__(f"{code_location}:: {msg}")


class SolarPanelType(enum.Enum):
    """
    Denotes the type of solar panel being modelled.

    - PV:
        Denotes a PV panel.

    - PV_T:
        Denotes a PV-T panel.

    - SOLAR_THERMAL
        Denotes a solar-thermal panel.

    """

    PV: str = "pv"
    PV_T: str = "pv_t"
    SOLAR_THERMAL: str = "solar_thermal"


def _conductive_air_heat_transfer_coefficient(wind_speed: float) -> float:
    """
    Calculates the conductive heat-transfer coefficient between PV panels and the air.

    The heat transfer coefficient is given by:
        h_air = 3.97(13) sWm^3K v_wind
              + 6.90(5) W/m^2K

    """

    return float(3.97 * wind_speed + 6.90)


def _sky_temperature(ambient_temperature: float) -> float:
    """
    Determines the radiative temperature of the sky.

    The "sky," as a black body, has a radiative temperature different to that of the
    surrounding air, or the ambient temperature. This function converts between them and
    outputs the sky's radiative temperature.

    Outputs:
        The radiative temperature of the "sky" in Kelvin.

    """

    return float(0.0552 * (ambient_temperature**1.5))


def _radiation_to_sky_coefficient(
    collector_emissivity: float, collector_temperature: float, sky_temperature: float
) -> float:
    """
    Returns the heat-transfer coefficient between a collector and the sky in W/m^2.

    Inputs:
        - collector_emissivity:
            The emissivity of the collector.
        - collector_temperature:
            The temperature of the collector, measured in Kelvin.
        - sky_temperature:
            The sky temperature, measured in Kelvin.

    Outputs:
        The heat-transfer coefficient between the collector and the sky, measured in
        Watts per meter squared.

    """
    return (
        STEFAN_BOLTZMAN_CONSTANT  # [W/m^2*K^4]
        * collector_emissivity
        * (collector_temperature**2 + sky_temperature**2)  # [K^2]
        * (collector_temperature + sky_temperature)  # [K]
    )


@dataclasses.dataclass
class PerformanceCurve:
    """
    Represents a performance curve for a solar-thermal collector.

    Solar-thermal collectors can be characterised by a performance curve,

        eta = eta_0 + c_1 * (T_c - T_a) / G + c_2 * (T_c - T_a)^2 / G,

    where `eta_0`, `c_1` and `c_2` give the zeroth-, first- and second-order
    coefficients which characterise the performance of the collector, `T_c` is the
    average temperature of the collector and `T_a` the ambient temperature, both
    measured in either degrees Kelvin or Celsius, but the same unit for each, and `G` is
    the solar irradiance, measured in Watts per meter squared.

    The attributes, `eta_0`, `c_1` and `c_2` are inherent properties of the collector
    and are contained within this class.

    .. attribute:: zeroth_order_cefficient
        The zeroth-order term for the performance curve.

    .. attribute:: first_order_cefficient
        The zeroth-order term for the performance curve.

    .. attribute:: second_order_cefficient
        The zeroth-order term for the performance curve.

    """

    zeroth_order_coefficient: float
    first_order_coefficient: float
    second_order_coefficient: float

    @property
    def eta_0(self) -> float:
        """
        Wrapper around the zeroth-order coefficient.

        Outputs:
            - The zeroth-order coefficient.

        """

        return self.zeroth_order_coefficient

    @property
    def c_1(self) -> float:
        """
        Wrapper around the first-order coefficient.

        Outputs:
            - The first-order coefficient.

        """

        return self.first_order_coefficient

    @property
    def c_2(self) -> float:
        """
        Wrapper around the second-order coefficient.

        Outputs:
            - The second-order coefficient.

        """

        return self.second_order_coefficient


@dataclasses.dataclass
class PVModuleCharacteristics:
    """
    Represents characteristcs of the PV module.

    .. attribute:: nominal_power
        Denoste the nominal power of the PV module in kWp.

    .. attribute:: reference_efficiency
        Denotes the reference efficiency of the PV module.

    .. attribute:: reference_temperature
        Denotes the reference temperature of the PV module, measured in degrees Kelvin.

    .. attribute:: thermal_coefficient
        Denotes the thermal coefficienct of the PV module.

    """

    nominal_power: float
    reference_efficiency: float
    _reference_temperature: float
    thermal_coefficient: float

    @property
    def reference_temperature(self) -> float:
        """
        Return the reference temperature in degrees Kelvin.

        Outputs:
            The reference temperature in degrees Kelvin.

        """

        return self._reference_temperature + ZERO_CELCIUS_OFFSET


class SolarPanel(abc.ABC, EmissableComponent):  # pylint: disable=too-few-public-methods
    """
    Represents a solar panel being considered.

    .. attribute:: area
        The area of the collector in meters squared, used to calculate the input power
        to the collector.

    .. attribute:: name
        The name of the panel being considered.

    .. attribite:: panel_type
        The type of panel being considered.

    """

    panel_type: SolarPanelType

    def __init__(
        self,
        area: float,
        cost: float,
        emissions: float,
        land_use: float,
        name: str,
        *,
        emissions_range: float = 0,
    ) -> None:
        """
        Instantiate a :class:`SolarPanel` instance.

        Inputs:
            - area:
                The surface area of the panel in meters squared.
            - cost:
                The cost of the solar panel per unit panel.
            - emissions:
                The emissions associated with the solar panel per unit panel.
            - land_use:
                The land occupied by the panel in meters squared.
            - name:
                The name to assign to the :class:`SolarPanel` in order to uniquely
                identify it.

        """

        self.area: float = area
        self.land_use: float = land_use

        super().__init__(cost, emissions, name, emissions_range=emissions_range)

    def __init_subclass__(cls, panel_type: SolarPanelType) -> None:
        """
        The init_subclass hook, run on instantiation of the :class:`SolarPanel`.

        Inputs:
            - panel_type:
                The type of panel being considered.

        Outputs:
            An instantiated :class:`SolarPanel` instance.

        """

        cls.panel_type = panel_type

        return super().__init_subclass__()

    @abc.abstractmethod
    def calculate_performance(
        self,
        ambient_temperature: float,
        logger: Logger,
        solar_irradiance: float,
        htf_heat_capacity: float | None,
        input_temperature: float | None,
        mass_flow_rate: float | None,
        wind_speed: float | None,
    ) -> Tuple[
        float | None,
        float | None,
        float | None,
        float | None,
    ]:
        """
        Abstract method for calculation of collector performance.

        Inputs:
            - ambient_temperature:
                The ambient temperature, measured in degrees Kelvin.
            - htf_heat_capacity:
                The heat capacity of the HTF entering the collector, measured in Joules
                per kilogram Kelvin (J/kgK).
            - input_temperature:
                The input temperature of the HTF entering the collector, measured in
                in degrees Kelvin.
            - logger:
                The :class:`logging.Logger` to use for the run.
            - mass_flow_rate:
                The mass-flow rate of HTF passing through the collector, measured in
                kilograms per second.
            - solar_irradiance:
                The solar irradiance incident on the surface of the collector, measured
                in Watts per meter squared.
            - wind_speed:
                The wind speed in meters per second.

        Outputs:
            - electrical_efficiency:
                The electrical efficiency of the PV panel.
            - output_temperature:
                The output temperature of the HTF leaving the collector.
            - reduced_temperature:
                The reduced temperature of the collector.
            - thermal_efficiency:
                The thermal efficiency of the collector.

        """

    @classmethod
    @abc.abstractmethod
    def from_dict(
        cls,
        logger: Logger,
        solar_inputs: dict[str, Any],
    ) -> Any:
        """
        Instantiate a :class:`SolarPanel` instance based on the input data.

        Inputs:
            - logger:
                The :class:`logging.Logger` to use for the run.
            - solar_inputs:
                The solar input data specific to this panel.

        """


class PVPanel(SolarPanel, panel_type=SolarPanelType.PV):
    """
    Represents a photovoltaic panel.

    .. attribute:: absorptivity
        The absorptivity of the panel, defined between 0 (none of the incident light is
        absorbed by the panel) to 1 (all of the incident light is absorbed by the
        panel).

    .. attribute:: emissivity
        The emissivity of the panel, defined between 0 (none of the incident light on
        the panel is emitted) to 1 (the panel radiates all incoming light).

    .. attribute:: pv_unit
        The unit of PV power being considered, defaulting to 1 kWp.

    .. attribute:: reference_efficiency
        The efficiency of the PV layer under standard test conditions.

    .. attribute:: reference_temperature
        The reference temperature of the PV layer of the panel, measured in degrees
        Kelvin.

    .. attribute:: thermal_coefficient
        The thermal coefficient of performance of the PV layer of the panel, measured in
        kelvin^(-1).

    """

    def __init__(
        self,
        absorptivity: float,
        area: float,
        cost: float,
        emissions: float,
        emissivity: float,
        land_use: float,
        name: str,
        pv_unit: float,
        reference_efficiency: float,
        reference_temperature: float,
        thermal_coefficient: float,
        *,
        emissions_range: float = 0,
    ) -> None:
        """
        Instantiate a :class:`PVPanel` instance.

        Inputs:
            - absorptivity:
                The absorptivity of the panel.
            - area:
                The surface area of the panel in meters squared.
            - cost:
                The cost of the :class:`PVPanel`.
            - emissions:
                The emissions associated with the :class:`PVPanel`.
            - emissivity:
                The emissivity of the panel.
            - land_use:
                The land occupied by the panel in meters squared.
            - name:
                The name to assign to the :class:`PVPanel` in order to uniquely identify
                it.
            - pv_unit:
                The output power, in Watts, of the PV layer of the panel per unit panel
                installed.
            - reference_efficiency:
                The reference efficiency of the panel, if required, otherwise `None`.
            - reference_temperature:
                The temperature, in degrees Celsius, at which the reference efficiency
                is defined, if required, otherwise `None`.
            - thermal_coefficient:
                The thermal coefficient of the PV layer of the panel, if required,
                otherwise `None`.
            - emissions_range:
                The range of carbon emissions associated with the component.

        """

        super().__init__(
            area,
            cost,
            emissions,
            land_use,
            name,
            emissions_range=emissions_range,
        )

        self.absorptivity = absorptivity
        self.emissivity = emissivity
        self.pv_unit: float = pv_unit
        self.reference_efficiency: float = reference_efficiency
        self._reference_temperature: float = reference_temperature
        self.thermal_coefficient: float = thermal_coefficient

    @property
    def reference_temperature(self) -> float:
        """
        Return the reference temperature in degrees Kelvin.

        Outputs:
            The reference temperature in degrees Kelvin.

        """

        return self._reference_temperature + ZERO_CELCIUS_OFFSET

    def __repr__(self) -> str:
        """
        Return a nice-looking representation of the panel.

        Outputs:
            - A nice-looking representation of the panel.

        """

        return (
            "PVPanel("
            + f"name={self.name}"
            + f", area={self.area:.2g}"
            + f", reference_efficiency={self.reference_efficiency:.2g}"
            + f", reference_temperature={self.reference_temperature:.2g} K"
            + f" ({self._reference_temperature:.2g} degC)"
            + f", thermal_coefficient={self.thermal_coefficient:.2g}"
            + ")"
        )

    @classmethod
    def from_dict(cls, logger: Logger, solar_inputs: dict[str, Any]) -> Any:
        """
        Instantiate a :class:`PVPanel` instance based on the input data.

        Inputs:
            - logger:
                The logger to use for the run.
            - solar_inputs:
                The solar input data for the panel.

        Outputs:
            A :class:`PVPanel` instance.

        """

        logger.debug("Attempting to create PVPanel from solar input data.")
        try:
            solar_inputs.pop("type")
        except KeyError:
            pass

        return cls(
            **solar_inputs,
        )

    def _average_module_temperature(
        self,
        ambient_temperature: float,
        logger: Logger,
        solar_irradiance: float,
        wind_speed: float,
    ) -> float:
        """
        Calculate the temperature of the PV module in Kelvin.

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

        Inputs:
            - ambient_temperature:
                The ambient temperature, measured in degrees Kelvin.
            - logger:
                The logger to use for the run.
            - solar_irradiance:
                The solar irradiance incident on the collector, measured in Watts per
                meter squared.
            - wind_speed:
                The wind speed in meters per second.

        Outputs:
            The average temperature of the PV module installed.

        """

        # Calculate variable values which remain constant throughout the iteration
        sky_temperature: float = _sky_temperature(ambient_temperature)

        # Setup inputs for the iterative loop
        best_guess_average_temperature: float = (
            ambient_temperature + (solar_irradiance / 800) * 20
        )
        solution_found: bool = False

        # Loop through until a solution is found
        logger.debug("Beginning iterative calculation for the PV module temperature.")
        _iter_var: int = 0
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
                    * (1 + self.thermal_coefficient * self.reference_temperature)
                )
                + _conductive_air_heat_transfer_coefficient(wind_speed)
                * ambient_temperature
                + radiation_to_sky_coefficient * sky_temperature
            ) / (
                radiation_to_sky_coefficient
                + _conductive_air_heat_transfer_coefficient(wind_speed)
                - self.absorptivity
                * self.thermal_coefficient
                * self.reference_efficiency
                * solar_irradiance
            )

            # Break if this average temperature is within the required precision.
            if (
                abs(best_guess_average_temperature - average_temperature)
                < TEMPERATURE_PRECISION
            ):
                solution_found = True
                logger.debug(
                    "Solution for PV-module temperature of %s found, returning.",
                    average_temperature,
                )

            best_guess_average_temperature = average_temperature
            _iter_var += 1
            logger.debug(
                "No solution found yet, re-iterating. Best guess: %s",
                best_guess_average_temperature,
            )

        return average_temperature

    def calculate_performance(
        self,
        ambient_temperature: float,
        logger: Logger,
        solar_irradiance: float,
        htf_heat_capacity: float | None = None,
        input_temperature: float | None = None,
        mass_flow_rate: float | None = None,
        wind_speed: float | None = None,
    ) -> Tuple[float | None, float | None, float | None, float | None]:
        """
        Calcuates the electrical performance of the collector.

        NOTE: The default temperature characteristics from
          Faiman, D. (2008). "Assessing the outdoor operating temperature of
          photovoltaic modules." Progress in Photovoltaics 16(4): 307-315.
        as implemented in the pvlib module are provided in a commented-out form should a
        developer or user which to use this calculation in-lieu of the in-house
        calculation.

        Should non-Silicon panels wish to be considered, the default valuesfor the pvlib
        calculation should be altered accordingly or a more-accurate module-temperature
        module employed. For the in-house calculation, the parameters for the collector
        should be adjusted in the inputs file.

        Inputs:
            - ambient_temperature:
                The ambient temperature, measured in degrees Kelvin.
            - htf_heat_capacity:
                The heat capacity of the HTF entering the collector, measured in Joules
                per kilogram Kelvin (J/kgK).
            - input_temperature:
                The input temperature of the HTF entering the collector, measured in
                in degrees Kelvin.
            - logger:
                The :class:`logging.Logger` to use for the run.
            - mass_flow_rate:
                The mass-flow rate of HTF passing through the collector, measured in
                kilograms per second.
            - solar_irradiance:
                The solar irradiance incident on the surface of the collector, measured
                in Watts per meter squared.
            - wind_speed:
                The wind speed in meters per second.

        Outputs:
            - electrical_efficiency:
                The electrical efficiency of the PV panel.
            - `None`:
                There is no HTF outputted from a PV panel.
            - average_temperature:
                The average temperature of the collector.
            - `None`:
                There is no thermal efficiency associated with a PV panel.

        """

        # Determine the average temperature of the panel using empirical PVLib modules.
        # average_temperature = temperature.faiman(
        #     solar_irradiance, ambient_temperature - ZERO_CELCIUS_OFFSET
        # )

        # Determine the average temperature of the panel using an in-house heat-balance
        # calculation.
        if wind_speed is None:
            logger.error("Wind speed needed for PV temperature calculation.")
            raise ProgrammerJudgementFault(
                "PVPanel:calculate_perofrmance",
                "No wind speed provided for PV module average temperature calculation.",
            )

        average_temperature = self._average_module_temperature(
            ambient_temperature, logger, solar_irradiance, wind_speed
        )

        # Determine the fractional electrical performance and electrical efficiency of
        # the PV panel.
        electrical_efficiency = self.reference_efficiency * (
            1
            - self.thermal_coefficient
            * (ambient_temperature - self.reference_temperature)
        )

        return electrical_efficiency, None, average_temperature, None
