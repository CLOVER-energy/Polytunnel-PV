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

from math import radians

from ..pv_cell import get_irradiance, PVCell


class _BasePVCellTest(unittest.TestCase):
    """
    Base test for testing PV cells.

    """

    def setUp(self) -> None:
        super().setUp()

        self.pv_cell_azimuth = 180
        self.pv_cell_length = 5
        self.pv_cell_tilt = 25
        self.pv_cell_width = 0.5
        self.pv_cell_breakdown_voltage = 15
        self.pv_cell_reference_temperature = 25

        self.pv_cell = PVCell(
            self.pv_cell_azimuth,
            self.pv_cell_length,
            self.pv_cell_tilt,
            self.pv_cell_width,
            self.pv_cell_breakdown_voltage,
            self.pv_cell_reference_temperature,
        )


class TestPVCell(_BasePVCellTest):
    """
    Tests the functionality of the PVCell.

    """

    def setUp(self) -> None:
        """Setup mocks in common across the test cases."""

        super().setUp()

        self.pv_cell = PVCell(
            self.pv_cell_azimuth,
            self.pv_cell_length,
            self.pv_cell_tilt,
            self.pv_cell_width,
            self.pv_cell_breakdown_voltage,
            self.pv_cell_reference_temperature,
        )

    def test_eq(self) -> None:
        """Tests the equality of cells that have the same orientation."""

        # Check that two equal cells are equal
        other_pv_cell = PVCell(
            self.pv_cell_azimuth,
            self.pv_cell_length,
            self.pv_cell_tilt,
            self.pv_cell_width,
            self.pv_cell_breakdown_voltage,
            self.pv_cell_reference_temperature,
        )
        self.assertEqual(self.pv_cell, other_pv_cell)

        # Check that two cells with the same orientation but different parameters are
        # different
        same_tilt_pv_cell = PVCell(
            self.pv_cell_azimuth,
            self.pv_cell_length + 20,
            self.pv_cell_tilt,
            self.pv_cell_width + 20,
            self.pv_cell_breakdown_voltage + 20,
            self.pv_cell_reference_temperature + 20,
        )
        self.assertEqual(self.pv_cell, same_tilt_pv_cell)

        # Check that two cells with a different tilt are not equal
        different_tilt_pv_cell = PVCell(
            self.pv_cell_azimuth,
            self.pv_cell_length + 20,
            self.pv_cell_tilt + 20,
            self.pv_cell_width + 20,
            self.pv_cell_breakdown_voltage + 20,
            self.pv_cell_reference_temperature + 20,
        )
        self.assertNotEqual(self.pv_cell, different_tilt_pv_cell)

        # Check that two cells with a different azimuth are not equal
        different_azimuth_pv_cell = PVCell(
            self.pv_cell_azimuth + 20,
            self.pv_cell_length + 20,
            self.pv_cell_tilt,
            self.pv_cell_width + 20,
            self.pv_cell_breakdown_voltage + 20,
            self.pv_cell_reference_temperature + 20,
        )
        self.assertNotEqual(self.pv_cell, different_azimuth_pv_cell)

    def test_lt_cell_id(self) -> None:
        """Tests the equality of cells that have the same orientation."""

        # Check that two equal cells are equal
        other_pv_cell = PVCell(
            self.pv_cell_azimuth,
            self.pv_cell_length,
            self.pv_cell_tilt,
            self.pv_cell_width,
            self.pv_cell_breakdown_voltage,
            self.pv_cell_reference_temperature,
            _cell_id=self.pv_cell.cell_id + 1,
        )
        self.assertTrue(self.pv_cell < other_pv_cell)

    def test_lt_tilt(self) -> None:
        """Tests the less-than when based on tilt."""

        # Check that two equal cells are equal
        other_pv_cell = PVCell(
            self.pv_cell_azimuth,
            self.pv_cell_length,
            self.pv_cell_tilt + 1,
            self.pv_cell_width,
            self.pv_cell_breakdown_voltage,
            self.pv_cell_reference_temperature,
        )
        self.assertTrue(self.pv_cell < other_pv_cell)

    def test_azimuth_in_radians(self) -> None:
        """Check that the azimuth is correctly converted to radians."""

        self.assertEqual(radians(self.pv_cell_azimuth), self.pv_cell.azimuth_in_radians)

    def test_tilt_in_radians(self) -> None:
        """Check that the tilt is correctly converted to radians."""

        self.assertEqual(radians(self.pv_cell_tilt), self.pv_cell.tilt_in_radians)

    def test_bypass(self) -> None:
        """Test that the bypass function works."""

        self.pv_cell.bypass((first_bypass_voltage := -15))
        self.assertEqual(self.pv_cell.breakdown_voltage, first_bypass_voltage)

        self.pv_cell.bypass((second_bypass_voltage := -25))
        self.assertEqual(self.pv_cell.breakdown_voltage, second_bypass_voltage)


class TestGetIrradiance(_BasePVCellTest):
    """Tests the get-irradiance function on PV cells."""

    def test_when_dark(self) -> None:
        """
        Tests when the code should return zero due to the sun being below the horizon.

        """

        self.pv_cell = PVCell(
            self.pv_cell_azimuth,
            self.pv_cell_length,
            self.pv_cell_tilt,
            self.pv_cell_width,
            self.pv_cell_breakdown_voltage,
            self.pv_cell_reference_temperature,
        )

        # Test during an eclipse.
        self.assertEqual(
            get_irradiance(self.pv_cell, 100, 0, 180, 90),
            0,
            "`get_irradiance` should return 0 if during an eclipse",
        )

        # Test sun below the horizion.
        self.assertEqual(
            get_irradiance(self.pv_cell, 100, 1000, 180, 0),
            0,
            "`get_irradiance` should return 0 if the sun has set",
        )

    def test_no_dni(self) -> None:
        """Tests the case where DNI is not provided."""

        pv_cell = PVCell(
            self.pv_cell_azimuth,
            self.pv_cell_length,
            self.pv_cell_tilt,
            self.pv_cell_width,
            self.pv_cell_breakdown_voltage,
            self.pv_cell_reference_temperature,
        )

        self.assertEqual(
            get_irradiance(pv_cell, 100, 1000, 180, 45), 1303.0603598718658
        )

    def test_mainline(self) -> None:
        """Tests the mainline case where DNI is provided."""

        pv_cell = PVCell(
            self.pv_cell_azimuth,
            self.pv_cell_length,
            self.pv_cell_tilt,
            self.pv_cell_width,
            self.pv_cell_breakdown_voltage,
            self.pv_cell_reference_temperature,
        )

        self.assertEqual(
            get_irradiance(pv_cell, 100, 1000, 180, 45, 920), 971.544127095287
        )
