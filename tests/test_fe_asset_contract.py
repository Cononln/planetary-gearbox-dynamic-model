"""Small CSV-to-MAT contract checks for the exact ANSYS sensor assets."""

import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import unittest
from uuid import uuid4

import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src_py"))
from fe_ring_asset_io import validate_export_qa  # noqa: E402

RESULTS_ROOT = ROOT / "results"
FREQUENCIES_HZ = (1000.0, 1600.0)
RING_NODES = 84
SENSOR_NODES = 2
BASE_GLOBAL_DEG = 35.0
THREE_ANGLES = (0, 120, 240)
FOUR_ANGLES = (0, 90, 120, 240)


def write_csv(path, header, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def change_csv_cell(path, row_index, column_index, value):
    with path.open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.reader(stream))
    rows[row_index][column_index] = value
    write_csv(path, rows[0], rows[1:])


def write_export_qa(source, *, max120=0.5, nodes=SENSOR_NODES, modes=2):
    write_csv(source / "auto_sensor_check.csv", ("item", "value"), (
        ("ring_cx", 0.0), ("ring_cy", 0.0), ("n_s0", nodes),
        ("n_s120", nodes), ("n_s240", nodes), ("mean120", 0.25),
        ("max120", max120), ("mean240", 0.3), ("max240", 0.6),
        ("unit", "mm_to_m"),
    ))
    (source / "UNIFIED_FE_ASSET_EXPORT_FINISHED.txt").write_text(
        "OK: unified FE modal asset exported from the same RST result set.\n"
        f"NMODES={modes}, SENSOR_0={nodes}, SENSOR_120={nodes}, SENSOR_240={nodes}\n",
        encoding="utf-8",
    )


def write_source(source, angles, sensor_columns, swapped=False, numeric_labels=False):
    """Build two tiny modes; labels refer to clockwise offsets from 35 degrees."""
    write_csv(
        source / "modal_summary.csv",
        ("mode_id", "frequency_Hz"),
        ((mode, frequency) for mode, frequency in enumerate(FREQUENCIES_HZ, 1)),
    )

    ring_with_frequency = sensor_columns == 10
    ring_header = (
        ("mode_id", "frequency_Hz", "node_id", "x_m", "y_m", "z_m", "ux_m", "uy_m", "uz_m")
        if ring_with_frequency else
        ("mode_id", "node_id", "x_m", "y_m", "z_m", "ux_m", "uy_m", "uz_m")
    )
    ring_rows = []
    for mode, frequency in enumerate(FREQUENCIES_HZ, 1):
        for node in range(RING_NODES):
            phi = 2.0 * math.pi * node / RING_NODES
            radial = (math.cos(phi), math.sin(phi))
            row = (mode, node + 1, 0.063 * radial[0], 0.063 * radial[1], 0.0,
                   mode * 1e-6 * radial[0], mode * 1e-6 * radial[1], 0.0)
            if ring_with_frequency:
                row = (mode, frequency, *row[1:])
            ring_rows.append(row)
    write_csv(source / "mode_shapes_ring_teeth.csv", ring_header, ring_rows)

    if sensor_columns == 9:
        sensor_header = (
            "mode_id", "sensor_id", "node_id", "x_m", "y_m", "z_m", "ux_m", "uy_m", "uz_m"
        )
    else:
        sensor_header = (
            "mode_id", "frequency_Hz", "angle_deg", "node_id", "x_m", "y_m",
            "z_m", "ux_m", "uy_m", "uz_m"
        )
    sensor_rows = []
    for mode, frequency in enumerate(FREQUENCIES_HZ, 1):
        for index, angle in enumerate(angles):
            actual_angle = {120: 240, 240: 120}.get(angle, angle) if swapped else angle
            phi = math.radians(BASE_GLOBAL_DEG - actual_angle)
            radial = (math.cos(phi), math.sin(phi))
            for local in range(SENSOR_NODES):
                radius = 0.05 + 0.0002 * local
                amplitude = expected_radial_shape(index, mode, local)
                values = (radius * radial[0], radius * radial[1], 0.001 * local,
                          amplitude * radial[0], amplitude * radial[1], 0.0)
                if sensor_columns == 9:
                    sensor_label = angle if numeric_labels else f"SENSOR_{angle}"
                    sensor_rows.append((mode, sensor_label, 1000 + index * 10 + local, *values))
                else:
                    sensor_rows.append((mode, frequency, angle, 1000 + index * 10 + local, *values))
    write_csv(source / "mode_shapes_sensor.csv", sensor_header, sensor_rows)


def expected_radial_shape(index, mode, local):
    return (10 * mode + 2 * index + local) * 1e-6


class FEAssetContractTest(unittest.TestCase):
    def setUp(self):
        RESULTS_ROOT.mkdir(exist_ok=True)
        self.fixture = RESULTS_ROOT / f"fe_asset_contract_{uuid4().hex}"
        self.fixture.mkdir()
        self.addCleanup(self.cleanup_fixture)

    def cleanup_fixture(self):
        if self.fixture.resolve().parent != RESULTS_ROOT.resolve():
            raise RuntimeError("Temporary FE fixture escaped the results directory")
        shutil.rmtree(self.fixture)

    def run_builder(self, source, output, count):
        script = (
            "prepare_ansys_exact_three_sensor_asset.py" if count == 3
            else "prepare_ansys_exact_four_sensor_asset.py"
        )
        return subprocess.run(
            [sys.executable, str(ROOT / "scripts" / script), "--source", str(source),
             "--output", str(output), "--ring-nodes", str(RING_NODES),
             "--sensor-nodes", str(SENSOR_NODES)],
            cwd=ROOT, capture_output=True, text=True, check=False,
        )

    def assert_asset(self, path, angles):
        asset = loadmat(path)
        sensor_count = len(angles)
        self.assertEqual(asset["frequency_Hz"].shape, (1, 2))
        self.assertEqual(asset["contact_angle_rad"].shape, (RING_NODES, 1))
        self.assertEqual(asset["contact_xyz_m"].shape, (RING_NODES, 3))
        self.assertEqual(asset["input_shape_normal"].shape, (RING_NODES, 2))
        self.assertEqual(asset["sensor_angle_local_deg"].shape, (1, sensor_count))
        self.assertEqual(asset["sensor_angle_global_deg"].shape, (1, sensor_count))
        self.assertEqual(asset["sensor_shape_radial"].shape, (sensor_count, 2))
        self.assertEqual(asset["sensor_centre_xyz_m"].shape, (sensor_count, 3))
        np.testing.assert_allclose(asset["frequency_Hz"][0], FREQUENCIES_HZ)
        np.testing.assert_allclose(asset["sensor_angle_local_deg"][0], angles)
        np.testing.assert_allclose(asset["sensor_angle_global_deg"][0],
                                   [(BASE_GLOBAL_DEG - angle) % 360 for angle in angles], atol=1e-10)
        self.assertEqual(int(asset["sensor_mode_shapes_exact"][0, 0]), 1)
        self.assertEqual(int(asset["mass_normalization_si_verified"][0, 0]), 0)
        for index, angle in enumerate(angles):
            # Each declared clockwise label must match its own FE patch, not
            # merely the unordered set of three or four patch positions.
            centre = asset["sensor_centre_xyz_m"][index]
            centre_deg = math.degrees(math.atan2(centre[1], centre[0])) % 360
            expected_deg = (BASE_GLOBAL_DEG - angle) % 360
            angular_error = (centre_deg - expected_deg + 180) % 360 - 180
            self.assertAlmostEqual(angular_error, 0.0, places=9)
            for mode in (1, 2):
                expected = sum(expected_radial_shape(index, mode, local)
                               for local in range(SENSOR_NODES)) / SENSOR_NODES
                self.assertAlmostEqual(asset["sensor_shape_radial"][index, mode - 1],
                                       expected, places=12)

    def test_three_sensor_csv_formats(self):
        for columns in (9, 10):
            with self.subTest(columns=columns):
                source = self.fixture / str(columns) / "source"
                output = self.fixture / str(columns) / "output"
                source.parent.mkdir()
                source.mkdir()
                write_source(source, THREE_ANGLES, columns)
                result = self.run_builder(source, output, 3)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assert_asset(output / "ansys_ring_modal_transfer_asset_exact_three_sensor.mat",
                                  THREE_ANGLES)

    def test_swapped_three_sensor_labels_are_rejected(self):
        for columns in (9, 10):
            with self.subTest(columns=columns):
                source = self.fixture / str(columns) / "source"
                output = self.fixture / str(columns) / "output"
                source.parent.mkdir()
                source.mkdir()
                write_source(source, THREE_ANGLES, columns, swapped=True)
                result = self.run_builder(source, output, 3)
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("sensor", (result.stdout + result.stderr).lower())
                self.assertFalse((output / "ansys_ring_modal_transfer_asset_exact_three_sensor.mat").exists())

    def test_three_sensor_numeric_apdl_labels(self):
        source = self.fixture / "source"
        output = self.fixture / "output"
        source.mkdir()
        write_source(source, THREE_ANGLES, 9, numeric_labels=True)
        result = self.run_builder(source, output, 3)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assert_asset(output / "ansys_ring_modal_transfer_asset_exact_three_sensor.mat",
                          THREE_ANGLES)

    def test_four_sensor_csv_and_mat_shape(self):
        for columns in (9, 10):
            with self.subTest(columns=columns):
                source = self.fixture / str(columns) / "source"
                output = self.fixture / str(columns) / "output"
                source.parent.mkdir()
                source.mkdir()
                write_source(source, FOUR_ANGLES, columns)
                result = self.run_builder(source, output, 4)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assert_asset(output / "ansys_ring_modal_transfer_asset_exact_four_sensor.mat",
                                  FOUR_ANGLES)

    def test_swapped_four_sensor_labels_are_rejected(self):
        source = self.fixture / "source"
        output = self.fixture / "output"
        source.mkdir()
        write_source(source, FOUR_ANGLES, 10, swapped=True)
        result = self.run_builder(source, output, 4)
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("sensor", (result.stdout + result.stderr).lower())
        self.assertFalse((output / "ansys_ring_modal_transfer_asset_exact_four_sensor.mat").exists())

    def test_mode_id_sequence_and_pairing_are_checked_in_both_formats(self):
        cases = (
            (9, "modal_summary.csv", 2, 0, "3"),
            (9, "mode_shapes_ring_teeth.csv", 2, 0, "2"),
            (9, "mode_shapes_sensor.csv", 2, 0, "2"),
            (10, "mode_shapes_ring_teeth.csv", 2, 0, "2"),
            (10, "mode_shapes_sensor.csv", 2, 0, "2"),
            (10, "mode_shapes_sensor.csv", 2, 0, "1.5"),
        )
        for index, (columns, filename, row, column, value) in enumerate(cases):
            with self.subTest(columns=columns, file=filename, value=value):
                source = self.fixture / f"mode_{index}"
                source.mkdir()
                write_source(source, THREE_ANGLES, columns)
                change_csv_cell(source / filename, row, column, value)
                output = self.fixture / f"mode_output_{index}"
                result = self.run_builder(source, output, 3)
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("mode_id", result.stdout + result.stderr)
                self.assertFalse((output / "ansys_ring_modal_transfer_asset_exact_three_sensor.mat").exists())

    def test_embedded_frequency_must_match_modal_summary_and_pairing(self):
        for index, filename in enumerate(("mode_shapes_ring_teeth.csv", "mode_shapes_sensor.csv")):
            with self.subTest(file=filename):
                source = self.fixture / f"frequency_{index}"
                source.mkdir()
                write_source(source, THREE_ANGLES, 10)
                change_csv_cell(source / filename, 2, 1, "1600")
                output = self.fixture / f"frequency_output_{index}"
                result = self.run_builder(source, output, 3)
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("frequency", (result.stdout + result.stderr).lower())
                self.assertFalse((output / "ansys_ring_modal_transfer_asset_exact_three_sensor.mat").exists())

    def test_export_qa_pair_values_and_freshness(self):
        source = self.fixture / "qa"
        source.mkdir()
        write_source(source, THREE_ANGLES, 9)
        self.assertFalse(validate_export_qa(source)["available"])

        write_export_qa(source)
        report = validate_export_qa(
            source, expected_mode_count=2, expected_sensor_nodes=SENSOR_NODES,
            check_marker_freshness=True,
        )
        self.assertTrue(report["available"])
        self.assertAlmostEqual(report["max_mapping_error_mm"], 0.6)

        marker = source / "UNIFIED_FE_ASSET_EXPORT_FINISHED.txt"
        marker.unlink()
        with self.assertRaisesRegex(ValueError, "Incomplete FE export"):
            validate_export_qa(source)
        write_export_qa(source, max120=2.1)
        with self.assertRaisesRegex(ValueError, "2 mm limit"):
            validate_export_qa(source)
        write_export_qa(source, nodes=3)
        with self.assertRaisesRegex(ValueError, "expected 2"):
            validate_export_qa(source, expected_sensor_nodes=2)
        write_export_qa(source, modes=3)
        with self.assertRaisesRegex(ValueError, "expected 2"):
            validate_export_qa(source, expected_mode_count=2)
        write_export_qa(source)
        old_time = marker.stat().st_mtime - 10
        os.utime(marker, (old_time, old_time))
        with self.assertRaisesRegex(ValueError, "Stale FE export"):
            validate_export_qa(source, check_marker_freshness=True)
        self.assertTrue(validate_export_qa(source, check_marker_freshness=False)["available"])

    def test_three_sensor_builder_rejects_failed_apdl_mapping_qa(self):
        source = self.fixture / "bad_mapping"
        output = self.fixture / "bad_mapping_output"
        source.mkdir()
        write_source(source, THREE_ANGLES, 9)
        write_export_qa(source, max120=2.1)
        result = self.run_builder(source, output, 3)
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("2 mm limit", result.stdout + result.stderr)
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
