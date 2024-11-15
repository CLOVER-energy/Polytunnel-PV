#!/usr/bin/python3.10
########################################################################################
# pv_module.py - Curved-PV-module Python module.                                       #
#                                                                                      #
# Author: Yaar Safra, Ben Winchester                                                   #
# Copyright: Yaar Safra, 2023                                                          #
# Date created: 06/11/2023                                                             #
# License: Open source                                                                 #
# Time created: 10:34:00                                                               #
########################################################################################
"""
pv_module.py - The pv-module module for Polytunnel-PV.

This module provides functionality for the modelling of the PV module.

"""

import enum
import functools
import warnings

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from math import acos, asin, cos, degrees, isnan, radians, pi, sin
from multiprocessing import Pool
from numpy import matmul
from typing import Any, Callable, TypeVar, Type

from .bypass_diode import BypassDiode, BypassedCellString
from .pv_cell import CellType, PVCell

__all__ = (
    "CircularCurve",
    "CurvedPVModule",
    "CurveType",
    "ModuleType",
    "TYPE_TO_CURVE_MAPPING",
)

# BIFACIAL:
#   String used to parse information about cell bifaciality.
BIFACIAL: str = "Bifacial"

# Floating point precision:
#   The floating-point precision of the numbers to use when doing rotations.
FLOATING_POINT_PRECISION: int = 8

# Number of processes:
#   The number of processes to use for parallel computations.
NUMBER_OF_PROCESSES: int = 8


class ImplementationError(Exception):
    """Raised when an invalid implementation of the classes or methods is used."""


class UndergroundCellError(Exception):
    """Raised when a PV cell is underground."""


class CurveType(enum.Enum):
    """
    Denotes the type of curve. Useful in constructing curves.

    - CIRCULAR:
        Denotes a circular geometry.

    """

    CIRCULAR: str = "circular"


# Type variable for Curve and children.
_C = TypeVar(
    "_C",
    bound="Curve",
)

# TYPE_TO_CURVE_MAPPING:
#   Mapping between the curve type and curve instances.
TYPE_TO_CURVE_MAPPING: dict[CurveType, _C] = {}


@dataclass(kw_only=True)
class Curve(ABC):
    """
    Represents a curve around which the PV module curves.

    A curved surface, _e.g._, a polytunnel, has a given axis around which it curves (in
    the case of a polytunnel, this is its length) an equation for its curve (which may
    be a simple circle or a more complex shape like a parabola or hyperbola) and a
    length scale. These last two are wrapped up in a single function which takes a
    disaplcement along the curve and returns the angles at that point.

    Attributes:
        - curvature_azimuth:
            The azimuth angle for the curvature axis in degrees.
        - curvature_tilt:
            The tilt angle for the curvature axis in degrees.
        - get_angles_from_surface_disaplacement:
            A callable function which can return the azimuth and tilt at any point on
            the surface based on the distance from the central axis.

    """

    curvature_axis_azimuth: float = 180
    curvature_axis_tilt: float = 0
    name: str = ""
    _azimuth_rotation_matrix: list[list[float]] | None = None
    _tilt_rotation_matrix: list[list[float]] | None = None

    def __init_subclass__(cls, curve_type: CurveType) -> None:
        """
        Hook used to store the type of the curve.

        :param: **curve_type:**
            The type of the curve.

        """

        cls.curve_type = curve_type  # type: ignore [attr-defined]
        TYPE_TO_CURVE_MAPPING[curve_type] = cls

        super().__init_subclass__()

    @abstractmethod
    def get_angles_from_surface_displacement(
        self, displacement: float
    ) -> tuple[float, float]:
        """
        Abstract method that must be implemented in subclasses.
        Calculate the azimuth and zenith angles at a point along the curve.

        :param: **displacement:**
            The distance from the central axis.

        :returns:
            - A tuple, (azimuth, tilt), with angles in degrees.

        """

        raise NotImplementedError("This method must be implemented in subclasses")

    @property
    def zenith_rotation_angle(self) -> float:
        """
        The rotation about the zenith is not the azimuth this is calculated here.

        A user-specified azimuth angle gives and angle from North (y=0) which is the
        opposite to a rotation carried out.

        """

        return -self.curvature_axis_azimuth

    @property
    def azimuth_rotation_matrix(self) -> list[list[float]]:
        """
        The rotation matrix for an azimuth rotation.

        Returns:
            - A `list` of `list`s representing the matrix.

        """

        if self._azimuth_rotation_matrix is None:
            self._azimuth_rotation_matrix = [
                [
                    cos(radians(self.zenith_rotation_angle)),
                    -sin(radians(self.zenith_rotation_angle)),
                    0,
                ],
                [
                    sin(radians(self.zenith_rotation_angle)),
                    cos(radians(self.zenith_rotation_angle)),
                    0,
                ],
                [0, 0, 1],
            ]

        return self._azimuth_rotation_matrix

    @property
    def tilt_rotation_angle(self) -> float:
        """
        The rotation angle for the zenith rotation needs to be adjusted also.

        Rotations about the x axis result in a tilting in a southerly rather than a
        northerly direction. This is adjusted for here.

        """

        return -self.curvature_axis_tilt

    @property
    def tilt_rotation_matrix(self) -> list[list[float]]:
        """
        The rotation matrix for an azimuth rotation.

        Returns:
            - A `list` of `list`s representing the matrix.

        """

        if self._tilt_rotation_matrix is None:
            self._tilt_rotation_matrix = [
                [1, 0, 0],
                [
                    0,
                    cos(radians(self.tilt_rotation_angle)),
                    -sin(radians(self.tilt_rotation_angle)),
                ],
                [
                    0,
                    sin(radians(self.tilt_rotation_angle)),
                    cos(radians(self.tilt_rotation_angle)),
                ],
            ]

        return self._tilt_rotation_matrix

    def _get_rotated_angles_from_surface_normal(
        self, un_rotated_normal: list[float]
    ) -> tuple[float, float]:
        """
        Rotate the normal vector based on the orientation of the curve/Polytunnel.

        :param: **un_rotated_normal:**
            The surface normal vector in the un-rotated frame.

        :returns:
            - A `tuple` containing information about the new rotated normal vector:
                - The azimuth angle, in degrees,
                - THe tilt angle, in degrees.

        """

        # Rotate this normal vector based on the tilt and azimuth of the polytunnel.
        rotated_normal = matmul(
            self.azimuth_rotation_matrix,
            matmul(self.tilt_rotation_matrix, un_rotated_normal),
        )

        # Compute the new azimuth and tilt angles based on these rotations.
        # The tilt angle is simply the z component of the vector.
        tilt_angle = round(acos(rotated_normal[2]), FLOATING_POINT_PRECISION)

        with warnings.catch_warnings():
            warnings.filterwarnings("error")
            try:
                x_y_plane_component = sin(tilt_angle)
            except (RuntimeError, RuntimeWarning, ZeroDivisionError):
                azimuth_angle = pi
            else:
                if round(rotated_normal[1], 6) == 0:
                    azimuth_angle = (
                        pi
                        if rotated_normal[0] == 0
                        else pi / 2 if (rotated_normal[0] > 0) else -pi / 2
                    )
                else:
                    arccos_angle = acos(
                        round(
                            rotated_normal[1] / x_y_plane_component,
                            FLOATING_POINT_PRECISION,
                        )
                    )
                    azimuth_angle = (
                        2 * pi - arccos_angle if rotated_normal[0] < 0 else arccos_angle
                    )

        # Check that the tilt is not out-of-bounds
        if tilt_angle > pi / 2:
            raise UndergroundCellError(
                f"A cell in the module has a tilt angle of {tilt_angle} radians, which "
                "is underground."
            )

        # Check that the azimuth angle is not `nan` and set it to South if so.
        if isnan(azimuth_angle):
            azimuth_angle = pi

        # Return these angles in degrees
        return degrees(azimuth_angle) % 360, degrees(tilt_angle)


@dataclass(kw_only=True)
class CircularCurve(Curve, curve_type=CurveType.CIRCULAR):
    """
    Represents a circular geometry. In this instance, a radius of curvature is required.

    Attributes:
        - radius_of_curvature:
            Represents the radius of curvature.

    """

    radius_of_curvature: float

    def get_angles_from_surface_displacement(
        self, displacement: float
    ) -> tuple[float, float]:
        """
        Calculate the azimuth and zenith angles at a point along the curve.

        :param: **displacement:**
            The distance from the central axis.

        :returns:
            - A tuple, (azimuth, tilt), with angles in degrees.

        """

        # Compute the zenith angle in radians based on the distance from the axis.
        zenith_angle: float = displacement / self.radius_of_curvature

        # Don't both calculating if the displacement is zero
        if displacement == 0:
            return (180, self.curvature_axis_tilt)

        # Compute the components of a unit normal vector with this zenith angle.
        un_rotated_normal: list[float] = [sin(zenith_angle), 0, cos(zenith_angle)]

        return self._get_rotated_angles_from_surface_normal(un_rotated_normal)


# Type variable for CurvedPVModule and children.
CPVM = TypeVar("CPVM", bound="CurvedPVModule")


class ModuleType(enum.Enum):
    """
    Denotes the type of the module.

    - THIN_FILM:
        Used to distinguish thin-film modules.

    """

    THIN_FILM: str = "thin_film"


@dataclass
class CurvedPVModule:
    """
    A curved PV module containing multiple solar cells.

    Attributes:
        - pv_cells:
            A `list` containing all the PV cells contained within the module.
        - offset_angle:
            The angle between the length of the PV module and the axis of the
            polytunnel.

    """

    pv_cells_and_cell_strings: list[BypassedCellString | PVCell]
    module_type: ModuleType
    name: str = ""
    offset_angle: float = 0

    def __post_init__(self) -> None:
        """
        Run post initialisation of the class constructor.

        Actions:
            - Check that the offset angle is within the allowed set of angles.
              NOTE: this will involve magic floats.

        """

        self._check_offset_angle_allowed()

    def __repr__(self) -> str:
        """The default representation of the module."""

        return (
            f"CurvedPVModule(name={self.name}, module_type={self.module_type.value}, "
            f"offset_angle={self.offset_angle:.3g}, "
            "pv_cells:\n\t{pv_cells}\n)".format(
                pv_cells="\n\t".join(
                    [str(cell) for cell in self.pv_cells_and_cell_strings]
                )
            )
        )

    def _check_offset_angle_allowed(self) -> None:
        """
        Check that the offset angle is allowed.

        Raises:
            - ImplementationErorr:
                Raised if the offset angle implementation chosen is not allowed.

        """

        if self.offset_angle not in (allowed_offset_angles := [0.0, 90.0]):
            raise ImplementationError(
                f"The offset angle of {self.offset_angle} degrees for the curved PV "
                "module is not allowed. Allowed values are "
                f"{', '.join(map(str, allowed_offset_angles))} degrees."
            )

    @property
    def pv_cells(self) -> list[PVCell]:
        """
        Return a `list` of _all_ the PV cells associated with the module.

        Returns:
            A `list` of all the PV cells, both bypassed and unbypassed.

        """

        pv_cells: list[PVCell] = []

        # Loop through the cell-or-string attribute and append or extend based on type
        for cell_or_bypassed_string in self.pv_cells_and_cell_strings:
            # If the entry is a pv cell, simply append.
            if isinstance(cell_or_bypassed_string, PVCell):
                pv_cells.append(cell_or_bypassed_string)
                continue

            # Otherwise, extend the list by the cells in the bypassed string.
            pv_cells.extend(cell_or_bypassed_string.pv_cells)

        return pv_cells

    @classmethod
    def thin_film_from_cell_number_and_dimensions(
        cls: Type[CPVM],
        cell_breakdown_voltage: float,
        cell_electrical_parameters: dict[str, float],
        cell_length: float,
        cell_spacing: float,
        cell_width: float,
        n_cells: int,
        *,
        bypass_diodes: list[BypassDiode],
        offset_angle: float,
        polytunnel_curve: Curve,
        module_centre_offset: float = 0,
        name: str = "",
    ) -> CPVM:
        """
        Instantiate a thin-film module based on the number of cells and dimensions.

        A thin-film module consists of a series of cells in a row:

            + | cell | cell | cell | cell | cell | cell | cell | cell | cell | -

        where each cell is a single strip.

        The cell dimensions and the number of cells determines the overall dimension of
        the module. The offset angle is the angle between the axis of the polytunnel and
        the **length** of the module. In this way, setting an offset angle of zero means
        that the cells are aligned along the axis of the polytunnel, whilst an offset
        angle of 90 (degrees) means that the cells are aligned perpendicular to the axis
        of the polytunnel.

        :param: **cell_breakdown_voltage:**
            The breakdown voltage of the PV cells, in Volts.

        :param: **cell_electrical_parameters:**
            Electrical parameters used to describe the IV curve of the cell.

        :param: **cell_length:**
            The length of the cells, _i.e._, the dimension parallel to the module
            length, given in meters.

        :param: **cell_spacing:**
            The space betweeh the cells, given in meters.

        :param: **cell_width:**
            The width of the cells, _i.e._, the dimension perpendicular to the
            module length, given in meters.

        :param: **n_cells:**
            The number of cells in the module.

        :param: **bypass_diodes:**
            The `list` of all bypass diodes present on the module.

        :param: **offest_angle:**
            The angle between the length of the module and the axis of the
            polytunnel.

        :param: **polytunnel_curve:**
            The curve defining the polytunnel.

        :param: **module_centre_offset:**
            The offset of the centre of the module from the centre of the curve, in
            meters.

        :returns:
            The instantiated `CurvedPVModule` instance.

        """

        def _cell_from_index(cell_index: int, cell_type: CellType) -> PVCell:
            """Construct a PV cell based on its index in the module."""

            # Compute the cell displacement based on its index and the module offset.
            cell_displacement = module_start_displacement + (
                (cell_index + 0.5) * cell_length + cell_index * cell_spacing
            ) * sin(radians(offset_angle))

            # Construct and return a PV cell.
            (
                cell_azimuth,
                cell_tilt,
            ) = polytunnel_curve.get_angles_from_surface_displacement(cell_displacement)
            return PVCell(
                azimuth=cell_azimuth,
                cell_type=cell_type,
                length=cell_length,
                tilt=cell_tilt,
                width=cell_width,
                breakdown_voltage=cell_breakdown_voltage,
                _cell_id=cell_index,
                **cell_electrical_parameters,  # type: ignore [arg-type]
            )

        # Determine whether the cells are mono- or bi-facial.
        cell_type: CellType = CellType(
            cell_electrical_parameters.pop(BIFACIAL, CellType.MONO_FACIAL.value)
        )

        # Determine the position of the start of the module.
        module_length = n_cells * cell_length + (n_cells - 1) * cell_spacing
        module_start_displacement = module_centre_offset - (module_length / 2) * sin(
            radians(offset_angle)
        )

        pv_cells: list[BypassedCellString | PVCell] = list(
            map(
                functools.partial(_cell_from_index, cell_type=cell_type),
                range(0, n_cells),
            )
        )

        # Bypass the cells accordingly:
        # Very first, confirm that none of the bypass diode ranges overlap.
        multiply_bypassed_cell_indicies = set()
        diode_ranges = [
            set(range(diode.start_index, diode.end_index)) for diode in bypass_diodes
        ]
        for diode_number, diode_range in enumerate(diode_ranges):
            for other_diode_range in diode_ranges[diode_number + 1 :]:
                multiply_bypassed_cell_indicies = multiply_bypassed_cell_indicies | (
                    diode_range & other_diode_range
                )

        if len(multiply_bypassed_cell_indicies) > 0:
            raise Exception(
                f"Bypass diodes for module {name} overlap. Multiply-bypassed cell "
                "indicies: "
                + ", ".join(
                    sorted([str(entry) for entry in multiply_bypassed_cell_indicies])
                )
            )

        # Go in reverse order, popping the cells into bypass diodes where appropriate.
        bypassed_cell_strings: list[BypassedCellString] = []
        for bypass_diode in reversed(
            sorted(bypass_diodes, key=lambda x: x.start_index)
        ):
            # Construct the bypassed cell string.
            bypassed_cell_strings.append(
                (
                    bypassed_cell_string := BypassedCellString(
                        bypass_diode=bypass_diode,
                        pv_cells=pv_cells[
                            bypass_diode.start_index : bypass_diode.end_index
                        ],
                    )
                )
            )

            # Remove these cells from the list.
            pv_cells = [
                cell for cell in pv_cells if cell not in bypassed_cell_string.pv_cells
            ]

        # Insert the bypassed cell strings into the list of PV cells and sort by cell
        # id.
        pv_cells.extend(bypassed_cell_strings)
        pv_cells = sorted(pv_cells, key=lambda cell: cell.cell_id)

        return cls(pv_cells, ModuleType.THIN_FILM, name, offset_angle)

    @classmethod
    def constructor_from_module_type(
        cls, module_type: ModuleType
    ) -> Callable[..., CPVM]:
        """
        Return the proper constructor based on the module type.

        :param: **module_type:
            The type of the PV module to use in the constructor.

        :returns: The appropriate constructor.

        """

        return {  # type: ignore [return-value]
            ModuleType.THIN_FILM: cls.thin_film_from_cell_number_and_dimensions
        }[module_type]


# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import pvlib as pv
# from pvlib import pvsystem
# from pvlib.irradiance import get_total_irradiance
# from scipy.constants import e as qe, k as kB
# import sys
# from datetime import datetime
# import seaborn as sns

# api_frame = pd.read_csv(
#     "/Users/yaarsafra/Library/CloudStorage/OneDrive-ImperialCollegeLondon/!Year 3 Imperial physics/3rd Year Project/Data/ninja_pv_43.6338_0.6031_uncorrected.csv",
#     skiprows=3,
#     delimiter=",",
# )
# api_frame["irradiance_direct"] = 1000 * api_frame["irradiance_direct"]
# api_frame["irradiance_diffuse"] = 1000 * api_frame["irradiance_diffuse"]

# """
# all functions below (above the SolarCell class) have been imported from
# the 'Calculating power loss from partial module shading' page in pvlib
# """


# def second_largest(list):
#     list.sort()
#     return list[-2]


# def simulate_full_curve(parameters, Geff, Tcell, ivcurve_pnts=1000, break_volt=-15):
#     """
#     Use De Soto and Bishop to simulate a full IV curve with both
#     forward and reverse bias regions.
#     """
#     # adjust the reference parameters according to the operating
#     # conditions using the De Soto model:
#     sde_args = pvsystem.calcparams_desoto(
#         Geff,
#         Tcell,
#         alpha_sc=parameters["alpha_sc"],
#         a_ref=parameters["a_ref"],
#         I_L_ref=parameters["I_L_ref"],
#         I_o_ref=parameters["I_o_ref"],
#         R_sh_ref=parameters["R_sh_ref"],
#         R_s=parameters["R_s"],
#     )
#     # sde_args has values:
#     # (photocurrent, saturation_current, resistance_series,
#     # resistance_shunt, nNsVth)

#     # Use Bishop's method to calculate points on the IV curve with V ranging
#     # from the reverse breakdown voltage to open circuit
#     kwargs = {
#         "breakdown_factor": parameters["breakdown_factor"],
#         "breakdown_exp": parameters["breakdown_exp"],
#         "breakdown_voltage": break_volt,  # THE MAXIMUM -VE VOLATAGE (WHERE THE VOLTAGE SPIKES UP)
#     }
#     v_oc = pv.singlediode.bishop88_v_from_i(
#         0.0, *sde_args, **kwargs
#     )  # CALCULATES THE VOLTAGE AT I = 0 (MAXIMUM +VE VOLTAGE OF IV CURVE)
#     # print('v_oc =', v_oc)
#     # print('breakdown voltage =', kwargs["breakdown_voltage"])
#     # ideally would use some intelligent log-spacing to concentrate points
#     # around the forward- and reverse-bias knees, but this is good enough:
#     vd = np.linspace(
#         0.99 * kwargs["breakdown_voltage"], v_oc, ivcurve_pnts
#     )  # GENERATES A LIST OF VOLTAGES BETWEEN THE MAXIMUM -VE AND +VE VOLTAGES
#     ivcurve_i, ivcurve_v, ivcurve_P = pv.singlediode.bishop88(
#         vd, *sde_args, **kwargs
#     )  # CALCULATES THE CURRENT AT THE VOLTAGE VALUES GIVEN

#     """CALCULATE THE IV CURVE BY CREATING A CURRENT LINSPACE INSTEAD TO MAKE PI CURVES EASIER TO ANALYSE"""
#     i_0, v_0, p_0 = pv.singlediode.bishop88(
#         0, *sde_args, **kwargs
#     )  # THE CURRENT, VOLTAGE AND POWER VALUES AS V = 0 FOR PI CURVE
#     i_values_PI = np.linspace(0, 10, ivcurve_pnts)
#     v_values_PI = pv.singlediode.bishop88_v_from_i(i_values_PI, *sde_args, **kwargs)
#     # print('MAX I =',i_0)
#     return pd.DataFrame(  # could also return the power at each point to the database, unsure if this makes things easier
#         {
#             "i_PI": i_values_PI,
#             "v_PI": v_values_PI,
#             "i": ivcurve_i,
#             "v": ivcurve_v,
#         }
#     )


# def plot_curves(dfs, labels):
#     """plot the forward- and reverse-bias portions of an IV curve"""
#     fig, axes = plt.subplots(1, 2, sharey=True, figsize=(15, 9))
#     index = 0
#     for df, label in zip(dfs, labels):
#         df.plot("v", "i", label=label, ax=axes[0])
#         df.plot("v", "i", label=label, ax=axes[1])
#         axes[0].set_xlim(right=0)
#         axes[0].set_ylim([0, 15])
#         axes[1].set_xlim([0, df["v"].max() * 1.5])
#         # print(df["i"][index])
#         # if df["i"]<0:
#         #     break
#         # i+=1
#     axes[0].set_ylabel("Current (A)", fontsize=40)
#     axes[0].set_xlabel("Voltage (V)", fontsize=40)
#     axes[1].set_xlabel("Voltage (V)", fontsize=40)
#     axes[0].legend_ = None
#     # axes[1].legend_ = None
#     axes[0].tick_params(axis="x", labelsize=30)
#     axes[0].tick_params(axis="y", labelsize=30)
#     axes[1].tick_params(axis="x", labelsize=30)
#     axes[1].tick_params(axis="y", labelsize=30)

#     # fig.suptitle(title)
#     fig.tight_layout()
#     axes[0].grid(alpha=0.5)
#     axes[1].grid(alpha=0.5)

#     axes[0].set_facecolor("gainsboro")
#     axes[1].set_facecolor("gainsboro")

#     # plt.rcParams.update({'font.size': 20})
#     # fig.sns.set_context("talk")
#     # plt.savefig("Figures/single cell iv curve", dpi=500)
#     return axes


# Vth = kB * (273.15 + 25) / qe


# SC = SolarCell
# # cell1 = SolarCell(x_width = 5, y_length = 1, cell_tilt=0, cell_azimuth = 90)
# # cell1.bypass()
# # cell1.get_iv_curve()


# class SolarModule:
#     def __init__(
#         self,
#         weather=0,
#         Rc=20,
#         N=3,
#         cell_length=3,
#         cell_width=1,
#         gap=0.5,
#         cell_azimuth=90,
#     ):
#         self._weather = weather  # adjust later (maybe split into different parameters)
#         self._Rc = Rc
#         self._N = int(N)
#         self._L = cell_length  # this length is not the arc length but the shorter, flat cell appriximation
#         self._W = cell_width
#         self._G = gap
#         self._cell_az = cell_azimuth

#     def get_Rc(self):
#         return self._Rc

#     def get_N(self):
#         return self._N

#     def get_length(self):
#         return self._L

#     def get_width(self):
#         return self._W

#     def get_gap(self):
#         return self._G

#     def get_curved_irradiance(
#         self, cell, solar_zenith, solar_azimuth, dni, dhi, curve=True
#     ):
#         # ghi = dhi + dni *np.cos(self._cell_tilt)
#         θs, φs, θc, φc = (
#             solar_zenith * (np.pi / 180),
#             solar_azimuth * (np.pi / 180),
#             cell._cell_tilt * (np.pi / 180),
#             cell._cell_azimuth * (np.pi / 180),
#         )
#         alignment_coef = (
#             np.sin(θc) * np.cos(φc) * np.sin(θs) * np.cos(φs)
#             + np.sin(θc) * np.sin(φc) * np.sin(θs) * np.sin(φs)
#             + np.cos(θc) * np.cos(θs)
#         )  # using the dot product of the vector normal to the surface of the cell and the vector of the sun relative to the cell
#         ghi = dhi + dni * alignment_coef
#         total_irrad = get_total_irradiance(
#             cell._cell_tilt,
#             cell._cell_azimuth,
#             solar_zenith,
#             solar_azimuth,
#             dni,
#             ghi,
#             dhi,
#         )  # (surface_tilt, surface_azimuth, solar_zenith, solar_azimuth, dni, ghi, dhi)
#         poa_global = total_irrad["poa_global"]
#         cell_angle = 2 * np.arcsin(self._L / (2 * self._Rc))
#         if solar_zenith > 90 or solar_zenith < -90:
#             poa_global = 0
#         if curve:
#             curvature_ratio = self._L / (
#                 self._Rc * cell_angle
#             )  # ratio between the length of the flat cell assumption to arclength of cell (<1)
#         else:
#             curvature_ratio = 1
#         return poa_global * curvature_ratio

#     def get_MPP(I, V):  # of a single cell (single IV curve)
#         P = [I * V for I, V in zip(I, V)]
#         MPP = np.max(P)
#         return MPP

#     def get_MPP_faster(I, V):  # of a single cell (single IV curve)
#         P = []
#         for i in range(len(I)):
#             power = I[i] * V[i]
#             P.append(power)
#             if power < 0:
#                 break

#         # P=[I[:10]*V[:10] for I[:10],V[:10] in zip(I[:10],V[:10])]
#         MPP = np.max(P)
#         return MPP

#     def get_iv_curve_module(self, solar_zen, solar_azi, cells, temperature):
#         curves = []
#         labels = []
#         for i in range(len(cells)):
#             curve = simulate_full_curve(
#                 cells[i].get_c_params(),
#                 SM.get_curved_irradiance(self, cells[i], solar_zen, solar_azi),
#                 temperature,
#                 break_volt=SC._cell_list[i].get_breakdown_voltage(),
#             )
#             # print(SM.get_curved_irradiance(self,cell[i]))
#             curves.append(curve)
#             # print(round(SC.get_cell_tilt(self[i]),2))
#             labels.append("angle =" + str(round(SC.get_cell_tilt(cells[i]), 2)))
#         # print(curves)
#         # print(labels)
#         axis = plot_curves(curves, labels)

#         # plt.rcParams.update({'font.size': 20})

#         return axis

#     def generate_flat_module(self, zenith=35, azimuth=180):
#         for i in range(self._N):
#             cell_i = SC(
#                 x_width=self._L,
#                 y_length=self._W,
#                 cell_tilt=zenith,
#                 cell_azimuth=azimuth,
#             )

#     def generate_curved_module(self):
#         cell_angle_dif = 2 * (
#             np.arcsin(self._L / (2 * self._Rc)) + np.arcsin(self._G / (2 * self._Rc))
#         )
#         cell_angle_dif *= 180 / np.pi  # convert to degrees
#         if self._N % 2 != 0:  # if odd number of cells
#             cell_i = SC(
#                 x_width=self._L,
#                 y_length=self._W,
#                 cell_tilt=0,
#                 cell_azimuth=self._cell_az,
#             )
#             for i in range((self._N - 1) // 2):
#                 # r =                                                           # TO DO: CALCULATE POSITION AND ADD TO CELL GENERATION LINE
#                 cell_angle = cell_angle_dif * (i + 1)
#                 if cell_angle > 90:
#                     print(
#                         "ERROR: one or more cells have an angle exceeds 90 degrees. The angle =",
#                         round(cell_angle, 2),
#                     )
#                     sys.exit()

#                 cell_i = SC(
#                     x_width=self._L,
#                     y_length=self._W,
#                     cell_tilt=cell_angle,
#                     cell_azimuth=self._cell_az,
#                 )
#                 cell_i = SC(
#                     x_width=self._L,
#                     y_length=self._W,
#                     cell_tilt=-cell_angle,
#                     cell_azimuth=self._cell_az,
#                 )
#                 # print('odd function, N =', self._N, ',i =', i, ',angle =', round(cell_angle,2))

#         else:
#             for i in range(self._N // 2):
#                 # r =                                                           # TO DO: CALCULATE POSITION AND ADD TO CELL GENERATION LINE
#                 cell_angle = cell_angle_dif / 2 + cell_angle_dif * i
#                 if cell_angle > 90:
#                     print(
#                         "ERROR: one or more cells have an angle exceeds 90 degrees. The angle =",
#                         round(cell_angle, 2),
#                     )
#                     sys.exit()
#                 cell_i = SC(
#                     x_width=self._L,
#                     y_length=self._W,
#                     cell_tilt=cell_angle,
#                     cell_azimuth=self._cell_az,
#                 )
#                 cell_i = SC(
#                     x_width=self._L,
#                     y_length=self._W,
#                     cell_tilt=-cell_angle,
#                     cell_azimuth=self._cell_az,
#                 )
#                 # print('even function, N =', self._N, ',i =', i, ',angle =', round(cell_angle,2))

#     def MPP_barchart(self, solar_zen=0, solar_azi=90, temperature=25):
#         cell = SC._cell_list
#         MPPs = []
#         tilts = []
#         for i in range(len(cell)):
#             curve = simulate_full_curve(
#                 cell[i]._c_params,
#                 SM.get_curved_irradiance(self, cell[i], solar_zen, solar_azi),
#                 temperature,
#                 break_volt=SC._cell_list[i].get_breakdown_voltage(),
#             )
#             MPP = SM.get_MPP(curve["i"], curve["v"])
#             MPPs.append(MPP)
#             tilts.append(round(SC.get_cell_tilt(cell[i]), 2))
#         plt.bar(tilts, MPPs, align="center", alpha=1, width=5)
#         # plt.xticks(y_pos, objects)
#         plt.ylabel("MPP (P)")
#         plt.xlabel("Angle (degrees)")

#         # plt.title("MPP Bar Chart of Module")
#         plt.savefig("Figures/7am MPP barchart", dpi=500)

#         plt.show()
#         # return axis

#     def irradiance_barchart(self, solar_zen=0, solar_azi=90):
#         cell = SC._cell_list
#         irradiances = []
#         tilts = []
#         for i in range(len(cell)):
#             irrad = SM.get_curved_irradiance(self, cell[i], solar_zen, solar_azi)
#             irradiances.append(irrad)
#             tilts.append(round(SC.get_cell_tilt(cell[i]), 2))

#         plt.bar(tilts, irradiances, align="center", alpha=1, width=5)
#         # plt.xticks(y_pos, objects)
#         plt.ylabel("Irradiance")
#         plt.title("Irradiance Bar Chart of Module")
#         plt.show()
#         # return axis

#     def power_current_graph(
#         self, cell=SC._cell_list, solar_zen=0, solar_azi=90, temperature=10
#     ):
#         # current = np.linspace(0,20,1000)
#         # curves=SC.get_iv_curve(SC._cell_list,show_axis=False) #use this line for obtaining multiple I V curves (implement the get_iv_curve_module method)
#         P_sum = []
#         all_P_values = []
#         negative_index = 1000
#         axis = plt.axes()
#         labels = []
#         for index in range(len(cell)):
#             curve = simulate_full_curve(
#                 cell[index]._c_params,
#                 SM.get_curved_irradiance(self, cell[index], solar_zen, solar_azi),
#                 temperature,
#                 break_volt=SC._cell_list[index].get_breakdown_voltage(),
#             )
#             P = []
#             I, V = curve["i_PI"], curve["v_PI"]
#             # print(V[51])
#             # print(I)
#             negative_check = True
#             for i in range(len(I)):
#                 power = I[i] * V[i]
#                 P.append(power)
#                 if index > 0:
#                     P_sum[
#                         i
#                     ] += power  # THE POWER VALUES ARE ADDED TO THE P_sum LIST FO RTHE REST OF THE CELLS
#                     # print('index =', index)
#                     # print('i =', i)
#                 if (
#                     power < -6 and negative_check == True
#                 ):  # THIS GETS RID OF THE ASYMPTOTES DISRUPTING THE FIGURE (ONLY PLOTS POINTS WITH A POWER > -6 WATTS)
#                     negative_index = i
#                     negative_check = False
#                     # print('NEGATIVE CHECKKK')
#             labels.append("angle =" + str(round(SC.get_cell_tilt(cell[index]), 2)))
#             if (
#                 index == 0
#             ):  # THE POWER VALUES OF THE FIRST CELL IS ADDED TO THE P_sum LIST
#                 P_sum += P
#                 # print('index =', index)
#             all_P_values += P
#             plt.plot(
#                 I[: negative_index + 1],
#                 P[: negative_index + 1],
#                 label=("angle =" + str(round(SC.get_cell_tilt(cell[index]), 2))),
#             )
#             # print(negative_index, P[negative_index-1])

#             plt.ylim(-2.5, 12)
#             plt.xlim(0, 10)

#             # plt.ylim(0, max_P*1.2)
#             plt.xlabel("Current (A)")
#             plt.ylabel("Power (W)")
#         # plt.savefig('Figures/single PI curve', dpi=500)
#         # print('length check: I =', len(I), ', V =', len(V), ', P =', len(P), ', P_sum =', len(P_sum))
#         # print('c =', c)
#         print("MPP =", max(P_sum))
#         max_P_sum_index = pd.Series(P_sum).idxmax()
#         max_I_sum = I[max_P_sum_index]
#         # max_P_sum_check = all_P_values[max_P_sum_index] + all_P_values[max_P_sum_index+1000] + all_P_values[max_P_sum_index+2000]
#         # print('max_P_sum_check =', max_P_sum_check)
#         # print(all_P_values[max_P_sum_index],all_P_values[max_P_sum_index])
#         y1_values = np.linspace(0, max(P_sum), 3)
#         x1_values = [max_I_sum, max_I_sum, max_I_sum]
#         x2_values = np.linspace(0, max_I_sum, 3)
#         y2_values = [max(P_sum), max(P_sum), max(P_sum)]
#         plt.plot(I, P_sum, color="black", ls="dotted")
#         # plt.plot(x1_values,y1_values,color='red',ls='dashed')
#         plt.plot(x2_values, y2_values, color="red", ls="dashed")
#         plt.plot([max_I_sum], [max(P_sum)], "gx", ms=20, color="red")
#         plt.gca().axhline(0, linestyle="--", color="grey")

#         plt.fill_between(I, P_sum, color="C3", alpha=0.3, hatch="//")
#         plt.rcParams.update({"font.size": 20})
#         plt.tight_layout()
#         # plt.subplots_adjust(bottom=0.8, right=0.9, top=0.9, left=0.8)
#         # plt.figure(figsize=(4, 3))

#         # plt.savefig('Figures/single PI curve2', dpi = 500)
#         return axis, max_I_sum, max(P_sum), P_sum

#     def MPP_BypDiodes_figure(
#         self, cells=SC._cell_list, solar_zen=0, solar_azi=90, temperature=10
#     ):
#         c = 0
#         if len(cells) % 2 == 1:
#             c = 1
#         cell_bypass_index = len(cells) - 1
#         max_P_sum_values = []
#         # max_I_sum_values =[]
#         number_of_bypassed = [0]
#         for index2 in range((len(cells) // 2) + c + 1):
#             # the last two cells in the _cell_list have the largest tilt
#             # print(cell_bypass_index)
#             if index2 != 0:
#                 cells[cell_bypass_index].bypass(-0.01)
#                 if cell_bypass_index != 0:
#                     cells[cell_bypass_index - 1].bypass(-0.01)
#                     number_of_bypassed.append((index2) * 2)
#                 else:
#                     number_of_bypassed.append((index2 * 2) - 1)
#                 cell_bypass_index -= 2
#             P_sum = []
#             P = []
#             current = []
#             for index in range(len(cells)):
#                 curve = simulate_full_curve(
#                     cells[index]._c_params,
#                     SM.get_curved_irradiance(
#                         self,
#                         cells[index],
#                         solar_zen,
#                         solar_azi,
#                         dni=200 * 5,
#                         ghi=220 * 5,
#                         dhi=20 * 5,
#                     ),
#                     temperature,
#                     break_volt=SC._cell_list[index].get_breakdown_voltage(),
#                 )  # PRODUCE IV VALUES
#                 # print(index2, cells[index])
#                 I, V = curve["i_PI"], curve["v_PI"]
#                 for i in range(len(I)):
#                     power = I[i] * V[i]
#                     P.append(power)
#                     if index > 0:
#                         P_sum[
#                             i
#                         ] += power  # THE POWER VALUES ARE ADDED TO THE P_sum LIST FOR THE REST OF THE CELLS
#                 if (
#                     index == 0
#                 ):  # THE POWER VALUES OF THE FIRST CELL IS ADDED TO THE P_sum LIST
#                     P_sum += P
#                 current.extend(I)
#                 # plt.savefig('Figures/single PI curve', dpi=500)
#                 # print('length check: I =', len(I), ', V =', len(V), ', P =', len(P), ', P_sum =', len(P_sum))
#             # max_P_sum_index = pd.Series(P_sum).idxmax()
#             # print(max(P_sum))
#             max_P_sum_values.append(max(P_sum))
#             # max_I_sum = I[max_P_sum_index]
#             # max_I_sum_values.append(max_I_sum)axis = plt.plot(I,P)
#             # plt.plot(current,P)
#             # plt.ylim(0, 2)
#             # plt.xlim(0, 2)
#             # plt.show()
#         print("MPP values:", max_P_sum_values)
#         plt.plot(number_of_bypassed, max_P_sum_values, "X", ls="dashed")
#         plt.ylim(
#             min(max_P_sum_values) * 0.999999999, max(max_P_sum_values) * 1.000000001
#         )
#         plt.xlabel("Number of Bypass Cells")
#         plt.ylabel("Power (W)")
#         plt.rcParams.update({"font.size": 20})
#         plt.tight_layout()
#         plt.savefig("Figures/Power vs bypass diodes 5am July", dpi=500)
#         plt.show()

#     def get_energy(
#         self,
#         solar_zen,
#         solar_azi,
#         temperature,
#         direct_ni,
#         diffuse_hi,
#         cells=SC._cell_list,
#         curved=True,
#     ):
#         irrad_module = 0
#         if curved:
#             # print('get_energy CURVED')
#             for i in range(len(cells)):
#                 irrad = SM.get_curved_irradiance(
#                     self,
#                     cells[i],
#                     solar_zen,
#                     solar_azi,
#                     dni=direct_ni,
#                     dhi=diffuse_hi,
#                     curve=curved,
#                 )
#                 irrad_module += irrad
#                 if (
#                     irrad < 0.00015
#                 ):  # to avoid runtime error and speed up whole operation
#                     sumP = [0, 0, 0]
#                     break
#                 else:
#                     data = simulate_full_curve(
#                         cells[i].get_c_params(),
#                         irrad,
#                         temperature,
#                         ivcurve_pnts=1000,
#                         break_volt=cells[i].get_breakdown_voltage(),
#                     )
#                     # (parameters, Geff, Tcell, ivcurve_pnts=1000,break_volt=-15)
#                     I, V = data["i_PI"], data["v_PI"]
#                     P = [I * V for I, V in zip(I, V)]
#                     if i == 0:
#                         sumP = P
#                     else:
#                         sumP = [sumP + P for sumP, P in zip(sumP, P)]
#         else:
#             # print('get_energy FLAT')
#             irrad = SM.get_curved_irradiance(
#                 self,
#                 cells[1],
#                 solar_zen,
#                 solar_azi,
#                 dni=direct_ni,
#                 dhi=diffuse_hi,
#                 curve=curved,
#             )
#             irrad_module = irrad * self._N
#             # print('irrad =', irrad)
#             if irrad == 0:
#                 sumP = [0, 0, 0]
#             else:
#                 data = simulate_full_curve(
#                     cells[1].get_c_params(), irrad, temperature, ivcurve_pnts=1000
#                 )
#                 # (parameters, Geff, Tcell, ivcurve_pnts=1000,break_volt=-15)
#                 I, V = data["i_PI"], data["v_PI"]
#                 P = [I * V for I, V in zip(I, V)]
#                 sumP = np.array(P) * self._N
#         # plt.plot(I,sumP,'x')
#         # plt.ylim(0,15)
#         # plt.xlim(0,0.5)
#         index_max = pd.Series(sumP).idxmax()
#         # print(index_max)
#         E_module = (sumP[index_max] * 3600) / (3.6 * 10**6)  # convert J to kWh
#         return E_module, irrad_module

#     def Energy_daily_figure(
#         self,
#         api_frame=api_frame,
#         time_index=np.arange(0, len(api_frame), 1),
#         cells=SC._cell_list,
#         curve=True,
#     ):
#         # data=generate_dataframe(self, api_frame, time_index, lat=43.3005, long=1.0905)
#         daily_energy = []
#         daily_energy_sum = 0
#         daily_irrad_sum = 0
#         daily_irrad_list = []
#         for i in range(len(time_index)):
#             date_time = datetime.fromtimestamp(
#                 time_index[i] * 3600 + 1546304400, tz=None
#             )  # convert time_index to datetime
#             solar_angles = get_solar_angles(
#                 date_time.year,
#                 date_time.month,
#                 date_time.day,
#                 date_time.hour,
#                 date_time.minute,
#                 latitude=43.3005,
#                 longitude=1.0905,
#                 temp=api_frame["temperature"][time_index[i]],
#             )
#             energy, irrad_module = SM.get_energy(
#                 self,
#                 solar_angles[0],
#                 solar_angles[1],
#                 api_frame["temperature"][time_index[i]],
#                 api_frame["irradiance_direct"][time_index[i]],
#                 api_frame["irradiance_diffuse"][time_index[i]],
#                 cells=SC._cell_list,
#                 curved=curve,
#             )  # calculate the energy produced every hour by using get_energy() function
#             # append this energy to hourly_energy
#             if (i + 1) % 24 != 0:
#                 print(energy, irrad_module)
#                 daily_energy_sum += energy
#                 daily_irrad_sum += irrad_module
#             else:
#                 daily_energy_sum += energy
#                 daily_energy.append(daily_energy_sum)
#                 daily_irrad_sum += irrad_module
#                 daily_irrad_list.append(daily_irrad_sum)
#                 print(
#                     i,
#                     "daily energy sum =",
#                     daily_energy_sum,
#                     "irradiance sum:",
#                     daily_irrad_sum,
#                 )
#                 daily_energy_sum = 0
#                 daily_irrad_sum = 0

#             # print(i)
#         # data.insert(4, "Energy", hourly_energy, True) #i think this is not needed
#         print("length=", len(daily_energy))
#         ax = plt.plot(np.arange(0, len(daily_energy), 1), daily_energy, "x")
#         # plt.xlim()
#         # plt.ylim(min(daily_energy)*0.99,second_largest(daily_energy)*1.01)
#         plt.xlabel("Days of the year 2019")
#         plt.ylabel("Energy (kWh)")

#         # ax.set_xticks(np.arange(1,13))
#         # ticklabels=['a','b','c','d','a','b','c','d','a','b','c','d']
#         # ax.set_xticklabels(ticklabels) #add monthlabels to the xaxis
#         plt.rcParams.update({"font.size": 15})

#         plt.tight_layout()
#         # plt.savefig('Figures/Power vs bypass diodes 5am July', dpi=500)
#         plt.show()
#         return daily_energy
#         # sum the energy for every day and add to


# sns.set_palette("PiYG", n_colors=2)
# sns.set_context("paper", font_scale=1.5)
# SM = SolarModule

# # E_day40=mod.Energy_daily_figure(time_index=np.arange(0*24,30*24,1))


# # E=SM.get_energy(mod,89,90,25,200,20,cells=SC._cell_list)
# # print(E)

# # cell1=SC(cell_tilt=0)
# # # cell2=SC(cell_tilt=20)
# # # print(cell2.get_c_params())
# # cell1.bypass(-10)
# # print(cell1)
# # print(cell1.get_c_params())
# # # print(cell2.get_c_params())
# # cell3=SC(cell_tilt=70)
# # cell3.bypass()


# # cell1.get_iv_curve()

# # mod.generate_curved_module()
# # SC._cell_list[0].bypass(breakdown_volt = -5)
# # mod.get_iv_curve_module(solar_zen=10, solar_azi=-90, cells=SC._cell_list, temperature=25)


# # c=0
# # if len(SC._cell_list)%2 == 1:
# #     c=1
# # cell_bypass_index=len(SC._cell_list)-1
# # max_P_sum_values =[]
# # # max_I_sum_values =[]
# # number_of_bypassed=[]
# # mod=SM(cell_length=1.1, N=10,Rc=5)
# # for index2 in range((len(SC._cell_list)//2)+c):
# #     #the last two cells in the _cell_list have the largest tilt
# #     SC._cell_list[cell_bypass_index].bypass(-15)
# #     if cell_bypass_index != 0:
# #         SC._cell_list[cell_bypass_index-1].bypass(-15)
# #         number_of_bypassed.append((index2+1)*2)
# #     else:
# #         number_of_bypassed.append((index2*2)+1)
# #     cell_bypass_index=-2
# #     mod.get_iv_curve_module(solar_zen=10, solar_azi=-90, cells=SC._cell_list, temperature=25)


# """single PI curve SLIDE 11"""
# # module=SM(N=1, Rc = 5)
# # SM.generate_curved_module(module)
# # module.power_current_graph(cell=SC._cell_list,solar_zen=0, solar_azi=90,temperature=25)     #currently only for one cell


# # change pallet: sns.set_palette("mako", n_colors=4)


# """Single IV curve SLIDE 10"""
# # cell1 = SolarCell(x_width = 5, y_length = 1, cell_tilt=0, cell_azimuth = 90)
# # print('cell list', SC._cell_list)

# # SC.get_iv_curve(SC._cell_list)


# # # # # cell2 = SolarCell(x_width = 3, y_length = 1, r_cell = np.array([2,0,0]), r_light = np.array([0,7,0]), cell_tilt=10, cell_azimuth = 0)
# # # # # cell3 = SolarCell(x_width = 3, y_length = 1, r_cell = np.array([2,0,0]), r_light = np.array([0,7,0]), cell_tilt=80, cell_azimuth = 0)
# # irrrad = SC.get_irradiance(cell1,dni=94/3.6, dhi=50/3.6,ghi = 98/3.6)

# # SM.get_curved_irradiance()


# """ MPP of 1 m^2 panel on an average day in London"""
# # cell1 = SolarCell(x_width = 5, y_length = 1, r_cell = np.array([-2,0,0]), cell_tilt=-70, cell_azimuth = 90)
# # irrrad = SC.get_irradiance(cell1,dni=94/3.6, dhi=50/3.6,ghi = 98/3.6)
# # curve = SC.get_iv_curve(SC._cell_list)
# # curve = simulate_full_curve(cell1._c_params, irrrad,11.6)
# # MPP = SM.get_MPP(curve['i'],curve['v'])
# # print(MPP)


# # cell_parameters = {
# #     "I_L_ref": 8.24,
# #     "I_o_ref": 2.36e-9,
# #     "a_ref": 1.3 * Vth,
# #     "R_sh_ref": 1000,
# #     "R_s": 0.00181,
# #     "alpha_sc": 0.0042,
# #     "breakdown_factor": 2e-3,
# #     "breakdown_exp": 3,
# #     "breakdown_voltage": -15,
# # }
# # kwargs = {
# #     "cell_parameters": cell_parameters,
# #     "poa_direct": 800,
# #     "poa_diffuse": 200,
# #     "Tcell": 25,
# # }
# # module_curve_full_sun = simulate_module(shaded_fraction=0, **kwargs)
# # module_curve_shaded = simulate_module(shaded_fraction=0.05, **kwargs)
# # ax = plot_curves([module_curve_full_sun, module_curve_shaded],
# #                  labels=['Full Sun', 'Shaded'],
# #                  title='Module-level reverse- and forward-biased IV curves')
# # plt.savefig('Full and shaded module combined', dpi = 200)

# """ figure of iv curves of multiple cells on a curved surface"""
# # print('initial no. of cells =', len(SC._cell_list))
# # module=SM(N=18)
# # SM.generate_curved_module(module)
# # SC.get_iv_curve(SC._cell_list)
# # print('final no. of cells =', len(SC._cell_list))


# # sns.set_palette("PiYG", n_colors=20)


# """ figure of MPP and irradiance bar chart of cells on a curved surface"""
# # mod=SM(cell_length=1.1, N=20,Rc=10)
# # mod.generate_curved_module()
# # mod.MPP_barchart(solar_zen=70, solar_azi=90)
# # mod.irradiance_barchart(solar_zen=, solar_azi=90)

# """ animation of MPP bar chart of cells on a curved surface from morning to evening on equator"""
# # mod=SM(cell_length=1.1, N=20,Rc=10)
# # mod.generate_curved_module()
# # for i in range(18):
# #     i-=9
# #     i=-i
# #     print(i)
# #     mod.MPP_barchart(solar_zen=i*10, solar_azi=90)


# """ importing the api (irradiances and temp) and the position of the Sun at a given
# time to get the MPPs of the different cells to create a cat (points) and strip plot"""


# def get_solar_angles(
#     year=2019,
#     month=1,
#     day=1,
#     hour=1,
#     minute=0,
#     latitude=43.3005,
#     longitude=1.0905,
#     temp=25,
# ):
#     date_time = "%g/%g/%g %g:%g" % (year, month, day, hour, minute)
#     date_time = datetime.strptime(
#         date_time, "%Y/%m/%d %H:%M"
#     )  # these two lines convert the date and time inputted into a striptime
#     date_time = date_time.timestamp()  # then convert to UTC timestamp
#     location = pv.location.Location(
#         latitude, longitude, tz="UTC", altitude=0, name=None
#     )  # set the location to a variable
#     solar_position = location.get_solarposition(
#         ((date_time) * 10**9),  # generate the angular postions of the Sun
#         pressure=None,
#         method="nrel_numpy",
#         temperature=temp,
#     )
#     # print(solar_position)
#     solar_azimuth = solar_position.azimuth[0]
#     solar_zenith = solar_position.apparent_zenith[
#         0
#     ]  # set the angular postions of the Sun to a variable
#     solar_angles = [solar_zenith, solar_azimuth]
#     return solar_angles


# def generate_times_index(
#     earliest, latest, step
# ):  # eariest and latest will be in a string format '%Y/%m/%s %H:%M', step will be in hours
#     selected_times_index = []
#     earliest = datetime.strptime(str(earliest), "%Y/%m/%d %H:%M").timestamp()
#     earl_index = (earliest - 1546304400) / 3600
#     # print('earliest index:', earl_index, 'ealiest timestamp:',earliest )
#     lastest = datetime.strptime(str(latest), "%Y/%m/%d %H:%M").timestamp()
#     late_index = (lastest - 1546304400) / 3600
#     # print('latest index:', late_index, 'latest timestamp:',lastest )
#     number_of_times = int((late_index - earl_index) / step)
#     # print('number of times:', number_of_times)
#     selected_times_index.append(earl_index)
#     for i in range(number_of_times + 1):
#         index = earl_index + i * step
#         selected_times_index.append(index)
#         if index > late_index:  # i dont think this is needed
#             print("ERROR: time index is greater that the latest time selected")
#             sys.exit()
#     return selected_times_index


# def generate_dataframe(
#     module, api_frame, time_index, lat=43.3005, long=1.0905, flat=False
# ):
#     dataframe = pd.DataFrame(
#         columns=["time index", "angle", "irradiance", "temperature", "MPP"]
#     )
#     # dataframe = pd.DataFrame(columns=['Column 1', 'Column 2', 'Column 3'])
#     # dataframe.columns =['time index', 'angle', 'irradiance', 'temperature']
#     for i in range(len(time_index)):  # repeat for every time index selected
#         date_time = datetime.fromtimestamp(
#             time_index[i] * 3600 + 1546304400, tz=None
#         )  # convert time_index to datetime
#         solar_angles = get_solar_angles(
#             date_time.year,
#             date_time.month,
#             date_time.day,
#             date_time.hour,
#             date_time.minute,
#             lat,
#             long,
#             api_frame["temperature"][time_index[i]],
#         )  # get solar zenith and azimuth
#         for j in range(
#             module.get_N()
#         ):  # repeat for every cell in the module to find the irradiance and add to the data frame
#             irrad = SM.get_curved_irradiance(
#                 module,
#                 SC._cell_list[j],
#                 solar_angles[0],
#                 solar_angles[1],
#                 dni=api_frame["irradiance_direct"][time_index[i]],
#                 dhi=api_frame["irradiance_diffuse"][time_index[i]],
#             )  # the dni and the dhi values are all zero!!!!
#             # print(irrad)
#             curve = simulate_full_curve(
#                 SC._cell_list[j].get_c_params(),
#                 irrad,
#                 api_frame["temperature"][time_index[i]],
#             )
#             # print(curve["v_PI"])
#             # plt.plot(curve["i_PI"][:10],curve["v_PI"][:10] ,'x')
#             # plt.ylim(-1,2)
#             # plt.xlim(0,1)
#             # print(max(curve["v_PI"]))
#             MPP = SM.get_MPP_faster(curve["i_PI"], curve["v_PI"])
#             print(
#                 "i=",
#                 i,
#                 "j=",
#                 j,
#                 "irrad=",
#                 irrad,
#                 "dni=",
#                 api_frame["irradiance_direct"][time_index[i]],
#                 "dhi=",
#                 api_frame["irradiance_diffuse"][time_index[i]],
#                 "MPP=",
#                 MPP,
#             )
#             new_row = {
#                 "time index": time_index[i],
#                 "angle": SC._cell_list[j].get_cell_tilt(),
#                 "irradiance": irrad,
#                 "temperature": api_frame["temperature"][time_index[i]],
#                 "MPP": MPP,
#             }  # create a row to add to the dataframe of the data calculated
#             dataframe.loc[len(dataframe)] = new_row  # add the new row to the dataframe
#     return dataframe


# """----IV curves of module-----"""
# # date_time = datetime.fromtimestamp(
# #     5 * 3600 + 182 * 24 * 3600 + 1546304400, tz=None
# # )  # convert time_index to datetime
# # solar_angles = get_solar_angles(
# #     date_time.year,
# #     date_time.month,
# #     date_time.day,
# #     date_time.hour,
# #     date_time.minute,
# #     latitude=43.3005,
# #     longitude=1.0905,
# #     temp=25
# # )  # get solar zenith and azimuth
# # SC._cell_list = sorted(SC._cell_list)
# # mod=SM(cell_length=1.1, N=20,Rc=10)


# # SM.generate_curved_module(mod)
# # # for i in range(len(SC._cell_list)):
# # #     SC._cell_list[i].bypass(-1)
# # axis = mod.get_iv_curve_module(cells = sorted(SC._cell_list), solar_zen = solar_angles[0],solar_azi=solar_angles[1], temperature=25)

# # plt.savefig('Figures/iv curves_6am_N20_0107', dpi = 500)

# """"----PI curves----"""
# # date_time = datetime.fromtimestamp(
# #     5 * 3600 + 182 * 24 * 3600 + 1546304400, tz=None
# # )  # convert time_index to datetime
# # solar_angles = get_solar_angles(
# #     date_time.year,
# #     date_time.month,
# #     date_time.day,
# #     date_time.hour,
# #     date_time.minute,
# #     latitude=43.3005,
# #     longitude=1.0905,
# #     temp=25
# # )  # get solar zenith and azimuth
# # SC._cell_list = sorted(SC._cell_list)
# # mod=SM(cell_length=1.1, N=18,Rc=10)
# # # module=SM(N=20, Rc = 10)
# # SM.generate_curved_module(mod)
# # # for i in range(len(SC._cell_list)):
# # #     SC._cell_list[i].bypass(-5)

# # axis, I_max, sum_P_max, P_sum=mod.power_current_graph(cell=sorted(SC._cell_list),solar_zen=solar_angles[0], solar_azi=solar_angles[1],temperature=25)     #currently only for one cell
# # plt.grid(alpha=0.5)
# # axis.set_facecolor('gainsboro')
# # plt.savefig('Figures/PI curves of module 6am 1st July', dpi=1500)
# # print(SC._cell_list)

# # mod.get_iv_curve_module(solar_zen=80, solar_azi=90, cells=SC._cell_list, temperature=10)
# # simulate_full_curve(SC._cell_list[0]._c_params,
# #     SM.get_curved_irradiance(mod, SC._cell_list[0], 80, 90),
# #     10,
# # )


# """!------------------------CAT PLOT-----------------------!"""
# # # import the irradiances and the temperatures at a given location and convert the time associated with them to a striptime and add this as a column to the data frame
# # api_frame = pd.read_csv(
# #     "/Users/yaarsafra/Library/CloudStorage/OneDrive-ImperialCollegeLondon/!Year 3 Imperial physics/3rd Year Project/Data/ninja_pv_43.6338_0.6031_uncorrected.csv",
# #     skiprows=3,
# #     delimiter=",",
# # )

# # UTC_time, local_time, electricity, dni, dhi, temp = (
# #     api_frame["time"],
# #     api_frame["local_time"],
# #     api_frame["electricity"],
# #     api_frame["irradiance_direct"],
# #     api_frame["irradiance_diffuse"],
# #     api_frame["temperature"],
# # )
# # time_stamps = []
# # for i in range(len(local_time)):
# #     time_stamp = datetime.strptime(local_time[i], "%d/%m/%Y %H:%M")
# #     time_stamps.append(time_stamp.timestamp())
# # api_frame["time_stamps"] = time_stamps


# # # generate the data frame for the cat plot using the methods created
# # mod = SM(cell_length=1, N=20, Rc=10)
# # mod.generate_curved_module()
# # time_indices = generate_times_index("2019/1/1 12:0", "2019/12/31 12:0", 24)
# # cat_frame = generate_dataframe(
# #     module=mod, api_frame=api_frame, time_index=time_indices, lat=43.3005, long=1.0905
# # )
# # cat_frame = cat_frame.round({"angle": 0})

# # print("-------...Generating catplot...-------")
# # sns.set_theme(
# #     rc={"figure.figsize": (20, 10)}
# # )  # make plot bigger to declutter the x axis lables
# # sns.set_context("paper", font_scale=1.5)

# # cat_facet_grid = sns.catplot(
# #     data=cat_frame,
# #     x="angle",
# #     y="MPP",
# #     hue="temperature",
# #     alpha=0.3,
# #     marker="D",
# #     zorder=1,
# #     palette="coolwarm",
# #     margin_titles="Cat plot irradiance of module (L=1,G=0.5,Rc=10,N=20) at 3pm everyday of the year 2019 in South of France",
# # )

# # # cat_facet_grid.fig.suptitle("Cat plot of irradiances (L=1,G=0.5,Rc=10,N=20) at 12pm")

# # cat_ax = plt.gca()
# # cat_ax.figure.set_size_inches(20, 10)
# # sns.set_context("paper", font_scale=3.5)

# # cat_facet_grid._legend.remove()
# # cat_ax.figure.colorbar(
# #     plt.cm.ScalarMappable(
# #         cmap="coolwarm",
# #         norm=plt.Normalize(
# #             cat_frame["temperature"].min(), cat_frame["temperature"].max(),

# #     )),
# #     ax=cat_ax,
# # )

# # # cat_ax.set(xlabel='Angles (degrees)', ylabel='Irradiance (W/m^2)', fontsize=10)
# # cat_ax.set_xlabel("Angle (degrees)",fontsize=30)
# # cat_ax.set_ylabel("Max Power Point (W)",fontsize=30)
# # cat_ax.set_ylim(-0.2,5)
# # # plt.grid(color='gray', alpha = 0.5)
# # # cat_ax.set_facecolor('gainsboro')
# # plt.savefig(f"Figures/Cat plot MPP of module", dpi=500)


# """"----MPP against no. of bypass diodes----"""


# # UTC_time, local_time, electricity, dni, dhi, temp = (
# #     api_frame["time"],
# #     api_frame["local_time"],
# #     api_frame["electricity"],
# #     api_frame["irradiance_direct"],
# #     api_frame["irradiance_diffuse"],
# #     api_frame["temperature"],
# # )

# # mod=SM(cell_length=1.1, N=20,Rc=10)
# # SM.generate_curved_module(mod)
# # hour = 14
# # day = 183
# # T = api_frame["temperature"][day*24+hour-1]
# # for i in range(1): #CREATING A LOOP DOES SOMETHING STRANGE TO THE MPP VALUES, THEY ARE ALL DIFFERENT IF INPUTED MANUALLY FOR THE SAME TIME AND THEY ARE CONSTANT W.R.T. THE NO. OF DIODES
# #     date_time = datetime.fromtimestamp(
# #     hour * 3600 + day * 24 * 3600 + 1546304400, tz=None
# # )  # convert time_index to datetime
# #     solar_angles = get_solar_angles(
# #     date_time.year,
# #     date_time.month,
# #     date_time.day,
# #     date_time.hour,
# #     date_time.minute,
# #     latitude=43.3005,
# #     longitude=1.0905,
# #     temp=T
# # )  # get solar zenith and azimuth
# # # SC._cell_list = sorted(SC._cell_list)
# # # module=SM(N=20, Rc = 10)
# #     mod.MPP_BypDiodes_figure(cells=SC._cell_list,solar_zen=solar_angles[0], solar_azi=solar_angles[1],temperature=T)     #currently only for one cell

# """"----Daily Energy throughout the year 2019 for curved module----"""

# # mod=SM(cell_length=1.1, N=20,Rc=10)
# # SM.generate_curved_module(mod)
# # E_day=mod.Energy_daily_figure()
# # E_day=[330822.62068216415, 342749.58170613315, 354841.5859087765, 354924.65171140205, 353694.8829027692, 359450.7918413245, 353497.15301727783, 333299.4161045042, 337401.22031562723, 353736.2522720855, 359309.20010676427, 343155.8707894677, 319197.839701285, 307128.3602765959, 347320.48129504937, 339907.2127095968, 336887.3683778086, 325266.8932071835, 340507.0038376632, 340833.6390355153, 344051.8548402089, 333133.71616233804, 335671.85729506495, 362330.35607369675, 326360.8737851864, 332125.7169767234, 336855.87488931837, 322021.4083908586, 316593.8884907338, 336179.5603846241, 320959.0932047963, 340602.9102775942, 331134.3479904841, 366508.2890373219, 346879.9274955578, 310497.28126441105, 336708.2871303104, 335806.4150814874, 344603.1887560526, 10301892886.836458, 336432.615570796, 347017.53180902934, 355852.4534288664, 382023.4255772442, 369673.9825555216, 371713.0153554602, 391484.9859080068, 398215.8393332682, 396667.7443358385, 402226.680627323, 399303.08931613533, 394497.18078468065, 391857.7672060429, 395384.55678941315, 398938.75080830674, 396701.13748807675, 393192.7340295858, 382585.8969394896, 402038.31047650427, 398608.32885621436, 399044.26263992384, 392488.58442184876, 401110.1453749447, 407820.8926007204, 376887.9785690354, 417548.2522287838, 412071.89985316695, 390078.24943445006, 402739.2442349872, 418876.473890676, 409030.88684363227, 419865.1378261203, 408313.59016421926, 395691.3377950723, 413717.179170948, 417187.4758834472, 427988.774495214, 423148.323670904, 432628.3892559932, 429069.83102525055, 417963.84861812554, 407219.06348316174, 408889.35058867483, 425668.6955732222, 460386.1335822417, 457190.30862194207, 473396.23467712634, 463973.06795735954, 467678.3221281874, 464354.85546611715, 457517.38660281437, 449652.47466290434, 469390.3975252649, 481752.53129395517, 467393.13030081696, 471994.4847635585, 479009.6353640479, 461934.5805885416, 477994.1717280346, 479600.9500883903, 468875.77920258226, 486368.2892886408, 488250.035385816, 471834.0255256646, 461607.8897288845, 477333.27194458403, 467657.57256989303, 458758.5496415468, 464601.5170448117, 475305.5633485456, 470441.0256654877, 467394.3130775849, 458563.7786970677, 476551.46131648886, 478830.6863504288, 483347.3408800733, 476412.15992690873, 487453.0918532776, 484375.7461770964, 488426.7928822089, 476371.42408890824, 483249.5558283268, 480649.8060338283, 488027.6036801268, 511306.3730203287, 500023.91816657083, 464674.29415595, 454889.97866109194, 507259.4431946264, 458627.77902495186, 480961.5939176134, 536547.4850608589, 539059.8940135193, 530889.8499416355, 531731.7400544435, 523086.35030627035, 506520.0380670879, 536502.5484638619, 530217.4741413867, 533730.1408410721, 535845.921859076, 526215.6629639397, 504362.09441437374, 494758.09730199736, 523546.84116467566, 514750.08993116085, 518641.9752409646, 534127.7399254831, 536254.544827566, 530477.7108236286, 515908.2694682, 501310.7819725607, 494084.93459468987, 504539.3900592371, 494619.6585954908, 527786.635650607, 533052.7437894174, 536265.8199933185, 535338.6962123401, 512359.3970800491, 522494.8829023009, 543411.164317806, 541374.8828928757, 521849.6224222506, 519561.61754649883, 522519.3958941324, 526872.3698304286, 505866.137365527, 482703.2657617453, 504099.24505002896, 499301.9572657847, 508697.7275133964, 519049.04920455226, 486396.1672624974, 489442.4970562453, 488699.4651122643, 469057.7371495701, 464513.5277211513, 489279.14366499917, 467665.6534886728, 488743.1867183189, 497552.4734083613, 491768.4317099985, 479574.6998777724, 473229.2280238428, 467725.8465379541, 480487.4105476561, 488742.6195741356, 470461.95448330266, 500660.24449028313, 502002.7628526637, 491583.3875556564, 492881.39969292516, 492289.89417653746, 495241.06257348374, 501815.27842901234, 496057.5382937383, 487552.22460676974, 503743.6041846381, 492168.96248029993, 481311.0139274655, 477743.27808727196, 460904.93986487534, 449057.0814004269, 454561.12722553604, 448394.4054356103, 464247.7765584678, 494825.9456055823, 505444.59476225794, 484640.70107689826, 497396.9030090915, 497603.41268307884, 495896.09796734527, 485520.1433646414, 456060.02258320217, 451227.59729332954, 458218.0576905223, 451443.0123650523, 441963.8975626404, 421728.2423192838, 412962.9734071615, 445204.4594090761, 453408.52299007995, 458039.3763222126, 458602.99509063334, 449505.6833016121, 451756.96380723093, 445369.72807438025, 434582.1124357327, 429264.1442018354, 444936.1149713666, 450008.63883489225, 448536.0838533299, 446906.2618933244, 430401.308531158, 420958.0113972268, 431592.02096960373, 419684.02454655175, 429614.0416734212, 431246.8694263699, 428098.2077985972, 417441.4117976693, 418697.4120202605, 427338.65148233087, 448147.4162058693, 442886.26349220274, 439004.751826639, 454444.89329926606, 454402.23177773634, 451310.00052290305, 449713.69383640005, 438893.73947997857, 444559.8423901332, 443187.0886432822, 433824.4562120823, 419998.5045429389, 415029.6104948813, 411164.97294604423, 392113.492632713, 388427.0846413113, 391485.7585782268, 402752.99990651646, 368548.6774926513, 366552.6181211739, 370506.97977809765, 387687.3600620515, 377983.7151132324, 391287.8626751876, 383763.2783999436, 384982.6589733007, 385471.49538510083, 371881.6927500637, 376000.6018186459, 375260.98315687437, 390048.0587711483, 389330.042716589, 382223.87881791365, 389229.5223052904, 379965.70606261876, 388311.5716747301, 372941.60459025064, 383387.72367426514, 389699.2115437355, 369060.4488625594, 348924.034611328, 357863.8563845226, 357569.88472640456, 383208.7051844959, 374799.80001888005, 368918.2231392619, 367146.9090287415, 369669.53842446744, 379504.02350417606, 392489.063064846, 316969.15676424635, 354795.1915670738, 376377.60080458096, 376574.2565106277, 366032.7527536025, 338757.402715859, 328114.08578003733, 341679.7001781, 343447.7863756298, 326221.1268861025, 310214.94683048385, 320790.94202310615, 336544.3204900124, 323597.6903011881, 320535.475271356, 331934.76657395833, 338250.015649083, 328580.7499888953, 330232.95804747264, 329865.1745481698, 329749.6980328684, 332640.5761043764, 334476.1323559745, 335591.5393304728, 330097.7308354531, 324247.5884769067, 314859.2878601249, 345768.4870283551, 352081.37449194165, 343870.9757857017, 338094.6223223498, 329073.30861701124, 315411.5272554186, 331890.3144005904, 322431.7271747897, 319208.1994249934, 323172.393845849, 315940.55001414183, 322738.14067043155, 312109.0042603083, 335779.1228600181, 336492.73311130615, 355954.33321709686, 331697.0804658759, 330600.64018588164, 333393.74007170927, 327507.45120295853, 310712.4346922036, 325451.2327364489, 330899.89574235585, 321055.3131721589, 296227.42607673915, 296989.7919948888, 313021.43734535744, 318228.1932340305, 305931.5552872913, 315939.35456444055, 319581.08525899204, 308170.4815210759, 318414.4027284329, 315464.5390831568, 315184.97490951454, 323182.9786964695, 308029.96809296444, 316489.38298673596, 318632.5386389429, 329829.83014249045, 333165.9679790495, 336291.4601311246, 334784.36067667615, 336291.70744316076]
# # sns.set_theme('paper')
# # ax=plt.figure()
# # plt.figure(figsize=(10,6))

# # plt.plot(np.arange(0,len(E_day),1),E_day,'x')
# # plt.ylim(0,0.3)
# # plt.xlabel('Days of the year 2019',fontsize=20)
# # plt.ylabel('Energy (kWh)',fontsize=20)
# # plt.grid(color='gainsboro')
# # # ax.set_xticks(np.arange(1,13))
# # # ticklabels=['a','b','c','d','a','b','c','d','a','b','c','d']
# # # ax.set_xticklabels(ticklabels) #add monthlabels to the xaxis
# # # plt.rcParams.update({'font.size': 15})
# # # ax.set_facecolor('white')

# # plt.tight_layout()
# # plt.savefig('Figures/Max Daily Energy curved panel with angle_coef', dpi=500, transparent=True)
# # plt.show()
# E_yearly_sum_curved = []
# for i in range(len(E_day)):
#     if E_day[i] < 0.3:
#         E_yearly_sum_curved.append(E_day[i])
# print(
#     "yearly curved energy output:",
#     np.sum(E_yearly_sum_curved),
#     "daily curved mean:",
#     np.mean(E_yearly_sum_curved),
# )


# """"----Daily Energy throughout the year 2019 for flat module----"""

# mod = SM(cell_length=1.1, N=11, Rc=11)
# SM.generate_flat_module(mod)
# # SM.get_curved_irradiance(mod, SC._cell_list, solar_zenith, solar_azimuth)
# E_day_flat = mod.Energy_daily_figure(curve=False)

# ax = plt.figure()
# plt.figure(figsize=(10, 6))

# plt.plot(np.arange(0, len(E_day_flat), 1), E_day_flat, "x")
# plt.ylim(0, 0.75)
# plt.xlabel("Days of the year 2019", fontsize=20)
# plt.ylabel("Energy (kWh)", fontsize=20)
# plt.grid(color="gainsboro")
# # ax.set_xticks(np.arange(1,13))
# # ticklabels=['a','b','c','d','a','b','c','d','a','b','c','d']
# # ax.set_xticklabels(ticklabels) #add monthlabels to the xaxis
# # plt.rcParams.update({'font.size': 15})
# # ax.set_facecolor('white')

# plt.tight_layout()
# plt.savefig(
#     "Figures/Max Daily Energy curved and flat panel same horizontal distance",
#     dpi=500,
#     transparent=True,
# )
# plt.show()

# E_yearly_sum_flat = []
# for i in range(len(E_day_flat)):
#     if E_day_flat[i] < 0.75:
#         E_yearly_sum_flat.append(E_day_flat[i])
# print(
#     "yearly flat energy output:",
#     np.sum(E_yearly_sum_flat),
#     "daily flat mean:",
#     np.mean(E_yearly_sum_flat),
# )
