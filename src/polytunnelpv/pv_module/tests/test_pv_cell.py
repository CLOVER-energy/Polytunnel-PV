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

from ..pv_cell import CellType, get_irradiance, PVCell


class _BasePVCellTest(unittest.TestCase):
    """
    Base test for testing PV cells.

    """

    def setUp(self) -> None:
        super().setUp()

        # Cell geometric parameters
        self.pv_cell_azimuth = 180
        self.pv_cell_length = 5
        self.pv_cell_tilt = 25
        self.pv_cell_type = CellType.MONO_FACIAL
        self.pv_cell_width = 0.5
        self.pv_cell_breakdown_voltage = 15
        self.pv_cell_reference_temperature = 25

        # Cell electrical parameters
        self.a_ref = 2
        self.reference_dark_current_density = 6.169936e-10
        self.reference_photocurrent_density = 9.791404
        self.reference_series_resistance = 0.255
        self.reference_shunt_resistance = 1800
        self.short_circuit_current_density_temperature_coefficient = -0.5

        self.pv_cell = PVCell(
            azimuth=self.pv_cell_azimuth,
            cell_type=self.pv_cell_type,
            length=self.pv_cell_length,
            tilt=self.pv_cell_tilt,
            width=self.pv_cell_width,
            breakdown_voltage=self.pv_cell_breakdown_voltage,
            reference_temperature=self.pv_cell_reference_temperature,
            a_ref=self.a_ref,
            reference_dark_current_density=self.reference_dark_current_density,
            reference_photocurrent_density=self.reference_photocurrent_density,
            reference_series_resistance=self.reference_series_resistance,
            reference_shunt_resistance=self.reference_shunt_resistance,
            short_circuit_current_density_temperature_coefficient=self.short_circuit_current_density_temperature_coefficient,
        )


class TestPVCellBasics(_BasePVCellTest):
    """
    Tests the functionality of the PVCell.

    """

    def setUp(self) -> None:
        """Setup mocks in common across the test cases."""

        super().setUp()

    def test_eq(self) -> None:
        """Tests the equality of cells that have the same orientation."""

        # Check that two equal cells are equal
        other_pv_cell = PVCell(
            azimuth=self.pv_cell_azimuth,
            cell_type=self.pv_cell_type,
            length=self.pv_cell_length,
            tilt=self.pv_cell_tilt,
            width=self.pv_cell_width,
            breakdown_voltage=self.pv_cell_breakdown_voltage,
            reference_temperature=self.pv_cell_reference_temperature,
            a_ref=self.a_ref,
            reference_dark_current_density=self.reference_dark_current_density,
            reference_photocurrent_density=self.reference_photocurrent_density,
            reference_series_resistance=self.reference_series_resistance,
            reference_shunt_resistance=self.reference_shunt_resistance,
            short_circuit_current_density_temperature_coefficient=self.short_circuit_current_density_temperature_coefficient,
        )
        self.assertEqual(self.pv_cell, other_pv_cell)

        # Check that two cells with the same orientation but different parameters are
        # different
        same_tilt_pv_cell = PVCell(
            azimuth=self.pv_cell_azimuth,
            cell_type=self.pv_cell_type,
            length=self.pv_cell_length + 20,
            tilt=self.pv_cell_tilt,
            width=self.pv_cell_width + 20,
            breakdown_voltage=self.pv_cell_breakdown_voltage + 20,
            reference_temperature=self.pv_cell_reference_temperature + 20,
            a_ref=self.a_ref,
            reference_dark_current_density=self.reference_dark_current_density,
            reference_photocurrent_density=self.reference_photocurrent_density,
            reference_series_resistance=self.reference_series_resistance,
            reference_shunt_resistance=self.reference_shunt_resistance,
            short_circuit_current_density_temperature_coefficient=self.short_circuit_current_density_temperature_coefficient,
        )
        self.assertEqual(self.pv_cell, same_tilt_pv_cell)

        # Check that two cells with a different tilt are not equal
        different_tilt_pv_cell = PVCell(
            azimuth=self.pv_cell_azimuth,
            cell_type=self.pv_cell_type,
            length=self.pv_cell_length + 20,
            tilt=self.pv_cell_tilt + 20,
            width=self.pv_cell_width + 20,
            breakdown_voltage=self.pv_cell_breakdown_voltage + 20,
            reference_temperature=self.pv_cell_reference_temperature + 20,
            a_ref=self.a_ref,
            reference_dark_current_density=self.reference_dark_current_density,
            reference_photocurrent_density=self.reference_photocurrent_density,
            reference_series_resistance=self.reference_series_resistance,
            reference_shunt_resistance=self.reference_shunt_resistance,
            short_circuit_current_density_temperature_coefficient=self.short_circuit_current_density_temperature_coefficient,
        )
        self.assertNotEqual(self.pv_cell, different_tilt_pv_cell)

        # Check that two cells with a different azimuth are not equal
        different_azimuth_pv_cell = PVCell(
            azimuth=self.pv_cell_azimuth + 20,
            cell_type=self.pv_cell_type,
            length=self.pv_cell_length + 20,
            tilt=self.pv_cell_tilt,
            width=self.pv_cell_width + 20,
            breakdown_voltage=self.pv_cell_breakdown_voltage + 20,
            reference_temperature=self.pv_cell_reference_temperature + 20,
            a_ref=self.a_ref,
            reference_dark_current_density=self.reference_dark_current_density,
            reference_photocurrent_density=self.reference_photocurrent_density,
            reference_series_resistance=self.reference_series_resistance,
            reference_shunt_resistance=self.reference_shunt_resistance,
            short_circuit_current_density_temperature_coefficient=self.short_circuit_current_density_temperature_coefficient,
        )
        self.assertNotEqual(self.pv_cell, different_azimuth_pv_cell)

    def test_lt_cell_id(self) -> None:
        """Tests the equality of cells that have the same orientation."""

        # Check that two equal cells are equal
        other_pv_cell = PVCell(
            azimuth=self.pv_cell_azimuth,
            cell_type=self.pv_cell_type,
            length=self.pv_cell_length,
            tilt=self.pv_cell_tilt,
            width=self.pv_cell_width,
            breakdown_voltage=self.pv_cell_breakdown_voltage,
            reference_temperature=self.pv_cell_reference_temperature,
            a_ref=self.a_ref,
            reference_dark_current_density=self.reference_dark_current_density,
            reference_photocurrent_density=self.reference_photocurrent_density,
            reference_series_resistance=self.reference_series_resistance,
            reference_shunt_resistance=self.reference_shunt_resistance,
            short_circuit_current_density_temperature_coefficient=self.short_circuit_current_density_temperature_coefficient,
            _cell_id=self.pv_cell.cell_id + 1,
        )
        self.assertTrue(self.pv_cell < other_pv_cell)

    def test_lt_tilt(self) -> None:
        """Tests the less-than when based on tilt."""

        # Check that two equal cells are equal
        other_pv_cell = PVCell(
            azimuth=self.pv_cell_azimuth,
            cell_type=self.pv_cell_type,
            length=self.pv_cell_length,
            tilt=self.pv_cell_tilt + 1,
            width=self.pv_cell_width,
            breakdown_voltage=self.pv_cell_breakdown_voltage,
            reference_temperature=self.pv_cell_reference_temperature,
            a_ref=self.a_ref,
            reference_dark_current_density=self.reference_dark_current_density,
            reference_photocurrent_density=self.reference_photocurrent_density,
            reference_series_resistance=self.reference_series_resistance,
            reference_shunt_resistance=self.reference_shunt_resistance,
            short_circuit_current_density_temperature_coefficient=self.short_circuit_current_density_temperature_coefficient,
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

        # Test during an eclipse.
        self.assertEqual(
            get_irradiance(self.pv_cell, 100, 0, 180, 0, 0),
            0,
            "`get_irradiance` should return 0 if during an eclipse",
        )

        # Test sun below the horizion.
        self.assertEqual(
            get_irradiance(self.pv_cell, 100, 1000, 180, 90, 0),
            0,
            "`get_irradiance` should return 0 if the sun has set",
        )

    @unittest.skip("DNI not implemented")
    def test_when_dark_no_dni(self) -> None:
        """
        Tests when the code should return zero due to the sun being below the horizon.

        """

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

    @unittest.skip("DNI not implemented")
    def test_no_dni(self) -> None:
        """Tests the case where DNI is not provided."""

        self.assertEqual(
            get_irradiance(self.pv_cell, 100, 1000, 180, 45), 1303.0603598718658
        )

    def test_mainline(self) -> None:
        """Tests the mainline case where DNI is provided."""

        self.assertEqual(
            get_irradiance(self.pv_cell, 100, 1000, 180, 45, 920), 971.544127095287
        )
