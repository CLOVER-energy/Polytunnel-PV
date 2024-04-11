#!/usr/bin/python3
########################################################################################
# test_pv_cell.py - Tests for the PV-cell module.                                      #
#                                                                                      #
# Author: Ben Winchester                                                               #
# Copyright: Ben Winchester, 2024                                                      #
# Date created: 11/04/2024                                                             #
# License: Open source                                                                 #
########################################################################################
"""
test_pv_cell.py - Tests for the PVCell code.

"""

import unittest

from math import degrees, pi
from unittest import mock

from ..pv_module import CircularCurve, CurvedPVModule


class _BaseCircularCurveTest(unittest.TestCase):
    """Test the circular-curve class."""

    def setUp(self) -> None:
        super().setUp()

        self.radius_of_curvature = 10
        self.tilted_north_south_circular_curve = CircularCurve(
            curvature_axis_azimuth=180,
            curvature_axis_tilt=10,
            radius_of_curvature=self.radius_of_curvature,
        )
        self.tilted_angled_circular_curve = CircularCurve(
            curvature_axis_azimuth=140,
            curvature_axis_tilt=10,
            radius_of_curvature=self.radius_of_curvature,
        )

    def _assert_angles_equal(
        self, first_set: tuple[float], second_set: tuple[float]
    ) -> None:
        """
        Used to enable `almostEqual` operation on `tuple` elements.

        Inputs:
            - first_set:
                The first set of angles to compare.
            - second_set:
                The second set of angles to compare.

        """

        self.assertAlmostEqual(first_set[0], second_set[0], places=5)
        self.assertAlmostEqual(first_set[1], second_set[1], places=5)

    def _tilt_angle_from_displacement(self, displacement: float) -> float:
        """
        Calculate the angle based on the displacement for this particular curve.

        Inputs:
            - displacement:
                The displacement, in meters.

        """

        return abs(degrees(displacement / self.curve.radius_of_curvature))


class TestSouthCurve(_BaseCircularCurveTest):
    """Tests a north-south curve."""

    def setUp(self) -> None:
        """Setup mocks in common across test cases."""

        super().setUp()

        self.curve = CircularCurve(
            curvature_axis_azimuth=180,
            curvature_axis_tilt=0,
            radius_of_curvature=self.radius_of_curvature,
        )

    def test_negatively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-3),
            (90, self._tilt_angle_from_displacement(-3)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-2),
            (90, self._tilt_angle_from_displacement(-2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-1),
            (90, self._tilt_angle_from_displacement(-1)),
        )

    def test_central_cell(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(0), (180, 0)
        )

    def test_positively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(1),
            (270, self._tilt_angle_from_displacement(1)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(2),
            (270, self._tilt_angle_from_displacement(2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(3),
            (270, self._tilt_angle_from_displacement(3)),
        )


class TestNorthCurve(_BaseCircularCurveTest):
    """Tests a north-south curve."""

    def setUp(self) -> None:
        """Setup mocks in common across test cases."""

        super().setUp()

        self.curve = CircularCurve(
            curvature_axis_azimuth=0,
            curvature_axis_tilt=0,
            radius_of_curvature=self.radius_of_curvature,
        )

    def test_negatively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-3),
            (270, self._tilt_angle_from_displacement(-3)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-2),
            (270, self._tilt_angle_from_displacement(-2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-1),
            (270, self._tilt_angle_from_displacement(-1)),
        )

    def test_central_cell(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(0), (180, 0)
        )

    def test_positively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(1),
            (90, self._tilt_angle_from_displacement(1)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(2),
            (90, self._tilt_angle_from_displacement(2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(3),
            (90, self._tilt_angle_from_displacement(3)),
        )


class TestEastFacingCurve(_BaseCircularCurveTest):
    """Tests a north-south curve."""

    def setUp(self) -> None:
        """Setup mocks in common across test cases."""

        super().setUp()

        self.curve = CircularCurve(
            curvature_axis_azimuth=90,
            curvature_axis_tilt=0,
            radius_of_curvature=self.radius_of_curvature,
        )

    def test_negatively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-3),
            (0, self._tilt_angle_from_displacement(-3)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-2),
            (0, self._tilt_angle_from_displacement(-2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-1),
            (0, self._tilt_angle_from_displacement(-1)),
        )

    def test_central_cell(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(0), (180, 0)
        )

    def test_positively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(1),
            (180, self._tilt_angle_from_displacement(1)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(2),
            (180, self._tilt_angle_from_displacement(2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(3),
            (180, self._tilt_angle_from_displacement(3)),
        )


class TestWestFacingCurve(_BaseCircularCurveTest):
    """Tests a north-south curve."""

    def setUp(self) -> None:
        """Setup mocks in common across test cases."""

        super().setUp()

        self.curve = CircularCurve(
            curvature_axis_azimuth=270,
            curvature_axis_tilt=0,
            radius_of_curvature=self.radius_of_curvature,
        )

    def test_negatively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-3),
            (180, self._tilt_angle_from_displacement(-3)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-2),
            (180, self._tilt_angle_from_displacement(-2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-1),
            (180, self._tilt_angle_from_displacement(-1)),
        )

    def test_central_cell(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(0), (180, 0)
        )

    def test_positively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(1),
            (0, self._tilt_angle_from_displacement(1)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(2),
            (0, self._tilt_angle_from_displacement(2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(3),
            (0, self._tilt_angle_from_displacement(3)),
        )


class TestNorthWestFacingCurve(_BaseCircularCurveTest):
    """Tests a north-south curve."""

    def setUp(self) -> None:
        """Setup mocks in common across test cases."""

        super().setUp()

        self.curve = CircularCurve(
            curvature_axis_azimuth=330,
            curvature_axis_tilt=0,
            radius_of_curvature=self.radius_of_curvature,
        )

    def test_negatively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-3),
            (240, self._tilt_angle_from_displacement(-3)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-2),
            (240, self._tilt_angle_from_displacement(-2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-1),
            (240, self._tilt_angle_from_displacement(-1)),
        )

    def test_central_cell(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(0), (180, 0)
        )

    def test_positively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(1),
            (60, self._tilt_angle_from_displacement(1)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(2),
            (60, self._tilt_angle_from_displacement(2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(3),
            (60, self._tilt_angle_from_displacement(3)),
        )


class TestNorthEastFacingCurve(_BaseCircularCurveTest):
    """Tests a north-south curve."""

    def setUp(self) -> None:
        """Setup mocks in common across test cases."""

        super().setUp()

        self.curve = CircularCurve(
            curvature_axis_azimuth=30,
            curvature_axis_tilt=0,
            radius_of_curvature=self.radius_of_curvature,
        )

    def test_negatively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-3),
            (300, self._tilt_angle_from_displacement(-3)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-2),
            (300, self._tilt_angle_from_displacement(-2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-1),
            (300, self._tilt_angle_from_displacement(-1)),
        )

    def test_central_cell(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(0), (180, 0)
        )

    def test_positively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(1),
            (120, self._tilt_angle_from_displacement(1)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(2),
            (120, self._tilt_angle_from_displacement(2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(3),
            (120, self._tilt_angle_from_displacement(3)),
        )


class TestSouthWestFacingCurve(_BaseCircularCurveTest):
    """Tests a north-south curve."""

    def setUp(self) -> None:
        """Setup mocks in common across test cases."""

        super().setUp()

        self.curve = CircularCurve(
            curvature_axis_azimuth=200,
            curvature_axis_tilt=0,
            radius_of_curvature=self.radius_of_curvature,
        )

    def test_negatively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-3),
            (110, self._tilt_angle_from_displacement(-3)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-2),
            (110, self._tilt_angle_from_displacement(-2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-1),
            (110, self._tilt_angle_from_displacement(-1)),
        )

    def test_central_cell(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(0), (180, 0)
        )

    def test_positively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(1),
            (290, self._tilt_angle_from_displacement(1)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(2),
            (290, self._tilt_angle_from_displacement(2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(3),
            (290, self._tilt_angle_from_displacement(3)),
        )


class TestSouthEastFacingCurve(_BaseCircularCurveTest):
    """Tests a north-south curve."""

    def setUp(self) -> None:
        """Setup mocks in common across test cases."""

        super().setUp()

        self.curve = CircularCurve(
            curvature_axis_azimuth=110,
            curvature_axis_tilt=0,
            radius_of_curvature=self.radius_of_curvature,
        )

    def test_negatively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-3),
            (20, self._tilt_angle_from_displacement(-3)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-2),
            (20, self._tilt_angle_from_displacement(-2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-1),
            (20, self._tilt_angle_from_displacement(-1)),
        )

    def test_central_cell(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(0), (180, 0)
        )

    def test_positively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(1),
            (200, self._tilt_angle_from_displacement(1)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(2),
            (200, self._tilt_angle_from_displacement(2)),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(3),
            (200, self._tilt_angle_from_displacement(3)),
        )
