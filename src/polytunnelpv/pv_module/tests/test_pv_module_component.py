#!/usr/bin/python3
########################################################################################
# test_pv_module_component.py - Tests for the PV-module component.                     #
#                                                                                      #
# Author: Ben Winchester                                                               #
# Copyright: Ben Winchester, 2024                                                      #
# Date created: 11/04/2024                                                             #
# License: Open source                                                                 #
########################################################################################
"""
test_pv_module_component.py - Tests for the component-level code.

"""

import math
import unittest

from ..pv_cell import relabel_cell_electrical_parameters
from ..pv_module import CircularCurve, CurvedPVModule, ImplementationError


class CurvedThinFilmPVModuleGeometry(unittest.TestCase):
    """Tests the curved PV-module code."""

    def setUp(self) -> None:
        """Setup mocks in common across the test cases."""

        # Electrical parameters that describe the cell.
        self._cell_electrical_parameters: dict[str, float | int | str] = (
            relabel_cell_electrical_parameters(
                {
                    # "Technology": "Mono-c-Si",
                    "Bifacial": 0,
                    # "STC": 360.214,
                    # "PTC": 324.9,
                    "A_c": 1.94,
                    # "Length": math.nan,
                    # "Width": math.nan,
                    "N_s": 72,
                    # "I_sc_ref": 9.79,
                    # "V_oc_ref": 47.2,
                    # "I_mp_ref": 9.26,
                    # "V_mp_ref": 38.9,
                    "alpha_sc": 0.005061,
                    "beta_oc": -0.163973,
                    # "T_NOCT": 46.3,
                    "a_ref": 2.009798,
                    "I_L_ref": 9.791404,
                    "I_o_ref": 6.169936e-10,
                    "R_s": 0.254823,
                    "R_sh_ref": 1777.348877,
                    # "Adjust": 8.897161,
                    "gamma_r": -0.4625,
                    # "BIPV": "N",
                    # "Version": "SAM 2018.11.11 r2",
                    # "Date": "1/3/2019",
                }
            )
        )
        self.curve = CircularCurve(
            curvature_axis_azimuth=180, curvature_axis_tilt=10, radius_of_curvature=10
        )

        super().setUp()

    def test_mainline(self) -> None:
        """Tests the mainline case."""

        # Test that the modules are arranged with increasing azimuthal orientation
        module = CurvedPVModule.thin_film_from_cell_number_and_dimensions(
            -15,
            self._cell_electrical_parameters,
            0.02,
            0.02,
            0.5,
            15,
            bypass_diodes=[],
            offset_angle=90,
            polytunnel_curve=self.curve,
            module_centre_offset=0,
        )

        # Check that 15 cells were made
        self.assertEqual(len(module.pv_cells), 15)

        # Check that the azimuthal angle increases for the cells
        for index, cell in enumerate(module.pv_cells[:-1]):
            self.assertTrue(module.pv_cells[index + 1].azimuth > cell.azimuth)

        # For a module along the axis, check that all of the cells have the same
        # orientation.
        module = CurvedPVModule.thin_film_from_cell_number_and_dimensions(
            -15,
            self._cell_electrical_parameters,
            0.02,
            0.02,
            0.5,
            15,
            bypass_diodes=[],
            offset_angle=0,
            polytunnel_curve=self.curve,
            module_centre_offset=0,
        )

        for cell in module.pv_cells[1:]:
            self.assertEqual(cell.azimuth, module.pv_cells[0].azimuth)

    def test_offset_angle_not_allowed(self) -> None:
        """Tests the case where the offset angle is not allowed."""

        with self.assertRaises(ImplementationError):
            # Test that the modules are arranged with increasing azimuthal orientation
            CurvedPVModule.thin_film_from_cell_number_and_dimensions(
                -15,
                self._cell_electrical_parameters,
                0.02,
                0.02,
                0.5,
                15,
                bypass_diodes=[],
                offset_angle=45,
                polytunnel_curve=self.curve,
                module_centre_offset=0,
            )
