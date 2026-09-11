"""问题 1 求解器的快速回归测试。"""

from __future__ import annotations

import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np

from src.q1_fvm import (
    DEFAULT_ENV_PATH,
    DEFAULT_PARAMETERS,
    DEFAULT_TEMPLATE_PATH,
    export_result_workbook,
    factor_tridiagonal,
    load_environment,
    make_grid,
    moisture_diffusivity,
    solve_heat,
    solve_mass,
)


class TestQ1Numerics(unittest.TestCase):
    def test_tridiagonal_solver_matches_dense_solver(self) -> None:
        lower = np.array([-1.0, -1.0, -1.0])
        diagonal = np.array([4.0, 4.0, 4.0, 4.0])
        upper = lower.copy()
        rhs = np.array([1.0, 2.0, 3.0, 4.0])
        dense = np.diag(diagonal) + np.diag(lower, -1) + np.diag(upper, 1)
        actual = factor_tridiagonal(lower, diagonal, upper).solve(rhs)
        np.testing.assert_allclose(actual, np.linalg.solve(dense, rhs), rtol=1e-13)

    def test_grid_volume_closes_cylinder_cross_section(self) -> None:
        _, radius, volume, _ = make_grid(80, DEFAULT_PARAMETERS.radius_m)
        self.assertAlmostEqual(radius[-1], DEFAULT_PARAMETERS.radius_m)
        self.assertAlmostEqual(
            float(volume.sum()), DEFAULT_PARAMETERS.radius_m**2 / 2.0, places=15
        )

    def test_diffusivity_is_positive_and_increases_with_moisture(self) -> None:
        values = moisture_diffusivity(np.array([0.1, 1.0, 2.55]))
        self.assertTrue(np.all(values > 0))
        self.assertTrue(np.all(np.diff(values) > 0))

    def test_closed_system_remains_at_initial_state(self) -> None:
        parameters = replace(
            DEFAULT_PARAMETERS,
            end_time_s=10.0,
            heat_transfer_w_m2_k=0.0,
            mass_transfer_m_s=0.0,
        )
        _, temperature, heat_diagnostics = solve_heat(
            20, 1.0, 10, lambda _time: 60.0, parameters=parameters
        )
        _, moisture, mass_diagnostics = solve_mass(
            20, 1.0, 10, lambda _time: 0.0, parameters=parameters
        )
        np.testing.assert_allclose(
            temperature, parameters.initial_temperature_c, atol=1e-12
        )
        np.testing.assert_allclose(
            moisture, parameters.initial_moisture_kg_kg, atol=1e-12
        )
        self.assertLess(heat_diagnostics["cumulative_relative_residual"], 1e-12)
        self.assertLess(mass_diagnostics["cumulative_relative_residual"], 1e-12)

    def test_official_environment_file_passes_validation(self) -> None:
        temperature, moisture, metadata = load_environment(DEFAULT_ENV_PATH)
        self.assertEqual(metadata["rows"], 241)
        self.assertAlmostEqual(temperature(0.0), 28.0)
        self.assertGreater(moisture(1800.0), 0.0)


class TestQ1Workbook(unittest.TestCase):
    def test_export_expands_official_template_and_rounds_results(self) -> None:
        times = np.array([1.0, 2.0, 3.0])
        radii = np.array([0.0, 0.1, 0.2])
        temperature = np.array(
            [
                [28.00001, 28.00006, 28.00009],
                [28.1, 28.2, 28.3],
                [28.4, 28.5, 28.6],
            ]
        )
        moisture = np.full((3, 3), 2.549987)
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_path = Path(temporary_directory) / "result1.xlsx"
            validation = export_result_workbook(
                DEFAULT_TEMPLATE_PATH,
                output_path,
                times,
                radii,
                temperature,
                moisture,
            )
            self.assertTrue(output_path.is_file())
            self.assertEqual(validation["expected_rows"], 4)
            self.assertEqual(validation["expected_columns"], 4)
            self.assertEqual(
                validation["sheets"]["温度"]["numeric_result_cells"], 9
            )


if __name__ == "__main__":
    unittest.main()
