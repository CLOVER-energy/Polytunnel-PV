#!/usr/bin/python3
########################################################################################
# test_pv_module.py - Tests for the PV-module module.                                  #
#                                                                                      #
# Author: Ben Winchester                                                               #
# Copyright: Ben Winchester, 2024                                                      #
# Date created: 11/04/2024                                                             #
# License: Open source                                                                 #
########################################################################################
"""
test_pv_module.py - Tests for the PVModule code.

"""

import unittest

from math import degrees

from ..pv_module import CircularCurve


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

        :param: **first_set:**
            The first set of angles to compare.

        :param: **second_set:**
            The second set of angles to compare.

        """

        self.assertAlmostEqual(first_set[0], second_set[0], places=3)
        self.assertAlmostEqual(first_set[1], second_set[1], places=3)

    def _tilt_angle_from_displacement(self, displacement: float) -> float:
        """
        Calculate the angle based on the displacement for this particular curve.

        :param: **displacement:**
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


class TestSlightlyTiltedSouthCurve(_BaseCircularCurveTest):
    """Tests a north-south curve."""

    def setUp(self) -> None:
        """Setup mocks in common across test cases."""

        super().setUp()

        self.curve = CircularCurve(
            curvature_axis_azimuth=180,
            curvature_axis_tilt=0.5,
            radius_of_curvature=self.radius_of_curvature,
        )

    def test_negatively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-3),
            (91.61589789067615, 17.195781234804805),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-2),
            (92.46499452691191, 11.469927509164924),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-1),
            (94.9707404963336, 5.751293051743691),
        )

    def test_central_cell(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(0), (180, 0.5)
        )

    def test_positively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(1),
            (265.02925950366637, 5.751293051743691),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(2),
            (267.53500547308806, 11.469927509164924),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(3),
            (268.38410210932386, 17.195781234804805),
        )


class TestTiltedSouthCurve(_BaseCircularCurveTest):
    """Tests a north-south curve."""

    def setUp(self) -> None:
        """Setup mocks in common across test cases."""

        super().setUp()

        self.curve = CircularCurve(
            curvature_axis_azimuth=180,
            curvature_axis_tilt=10,
            radius_of_curvature=self.radius_of_curvature,
        )

    def test_negatively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-3),
            (119.30805185402161, 19.809786583530162),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-2),
            (130.58436489497103, 15.16487503418409),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(-1),
            (149.98036440042233, 11.510607512619211),
        )

    def test_central_cell(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(0), (180, 10)
        )

    def test_positively_displaced_half(self) -> None:
        """Tests a flat, north-south circular curve."""

        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(1),
            (210.01963559957767, 11.510607512619211),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(2),
            (229.415635105029, 15.16487503418409),
        )
        self._assert_angles_equal(
            self.curve.get_angles_from_surface_displacement(3),
            (240.69194814597842, 19.809786583530162),
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
