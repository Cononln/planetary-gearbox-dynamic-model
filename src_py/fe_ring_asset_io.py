"""Read ANSYS modal CSV data for the three- and four-sensor asset builders."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import re

import numpy as np


SENSOR_IDS = ("SENSOR_0", "SENSOR_120", "SENSOR_240")
_FREQUENCY_RTOL = 1e-5  # Legacy CSVs may round embedded frequencies.


def _positive_mode_id(value: str, path: Path, row_index: int) -> int:
    try:
        number = float(value)
    except ValueError as exc:
        raise ValueError(f"Invalid mode_id in row {row_index} of {path}") from exc
    if not math.isfinite(number) or not number.is_integer() or number < 1:
        raise ValueError(f"Invalid mode_id {value!r} in row {row_index} of {path}")
    return int(number)


def _reference_frequencies(path: Path, mode_count: int) -> np.ndarray | None:
    summary = path.parent / "modal_summary.csv"
    if not summary.is_file():
        return None
    result = frequencies(summary)
    if result.size != mode_count:
        raise ValueError(
            f"{summary} contains {result.size} modes, expected {mode_count} for {path}"
        )
    return result


def _check_frequency(value: str, reference: float | None, path: Path, row_index: int) -> float:
    try:
        frequency = float(value)
    except ValueError as exc:
        raise ValueError(f"Invalid frequency in row {row_index} of {path}") from exc
    if not math.isfinite(frequency) or frequency <= 0:
        raise ValueError(f"Invalid frequency in row {row_index} of {path}")
    if reference is not None and not math.isclose(
        frequency, reference, rel_tol=_FREQUENCY_RTOL, abs_tol=1e-6
    ):
        raise ValueError(
            f"Frequency {frequency} in row {row_index} of {path} does not match "
            f"modal_summary ({reference})"
        )
    return frequency


def frequencies(path: Path) -> np.ndarray:
    """Return ordered positive modal frequencies in hertz."""
    values: list[float] = []
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        rows = csv.reader(f)
        next(rows, None)
        for row_index, row in enumerate(rows, start=2):
            if not row or not any(cell.strip() for cell in row):
                continue
            if len(row) < 2 or not row[1].strip():
                raise ValueError(f"Missing frequency in row {row_index} of {path}")
            mode_id = _positive_mode_id(row[0], path, row_index)
            if mode_id != len(values) + 1:
                raise ValueError(
                    f"modal_summary mode_id sequence in {path}: "
                    f"expected {len(values) + 1}, found {mode_id}"
                )
            values.append(_check_frequency(row[1], None, path, row_index))
    result = np.asarray(values)
    if result.size == 0 or not np.all(np.isfinite(result)) or np.any(result <= 0):
        raise ValueError(f"{path} does not contain positive modal frequencies")
    return result


def modal_blocks(path: Path, rows_per_mode: int, mode_count: int):
    """Read ring modal rows from either the legacy or frequency-inclusive CSV.

    The bundled APDL export writes ``mode,node,x,y,z,ux,uy,uz`` in SI units.
    The existing frequency-inclusive export writes
    ``mode,frequency,node,x,y,z,ux,uy,uz`` in N-mm-MPa. Values are converted
    to SI when the coordinate radius indicates millimetres.
    """
    xyz = np.empty((rows_per_mode, 3))
    displacement = np.empty((mode_count, rows_per_mode, 3))
    count = 0
    frequency_inclusive = None
    reference = _reference_frequencies(path, mode_count)
    embedded_frequency = np.full(mode_count, np.nan)
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        rows = csv.reader(f)
        next(rows, None)
        for row_index, row in enumerate(rows, start=2):
            if not row or not row[0].strip():
                continue
            if len(row) not in (8, 9):
                raise ValueError(f"Expected 8 or 9 fields in modal row of {path}")
            if frequency_inclusive is None:
                frequency_inclusive = len(row) == 9
            elif frequency_inclusive != (len(row) == 9):
                raise ValueError(f"Mixed modal row formats in {path}")
            mode, local = divmod(count, rows_per_mode)
            if mode >= mode_count:
                raise ValueError(f"Unexpected extra row in {path}")
            mode_id = _positive_mode_id(row[0], path, row_index)
            if mode_id != mode + 1:
                raise ValueError(
                    f"Ring mode_id block sequence in {path}: expected {mode + 1}, "
                    f"found {mode_id} in row {row_index}"
                )
            if frequency_inclusive:
                expected_frequency = (
                    reference[mode] if reference is not None else
                    embedded_frequency[mode] if local else None
                )
                embedded_frequency[mode] = _check_frequency(
                    row[1], expected_frequency, path, row_index
                )
                point = np.asarray(row[3:6], dtype=float)
                values = np.asarray(row[6:9], dtype=float)
            else:
                point = np.asarray(row[2:5], dtype=float)
                values = np.asarray(row[5:8], dtype=float)
            if mode == 0:
                xyz[local] = point
            elif not np.allclose(point, xyz[local], rtol=0.0, atol=1e-12):
                raise ValueError(f"Node ordering changed in {path}")
            displacement[mode, local] = values
            count += 1
    expected = mode_count * rows_per_mode
    if count != expected:
        raise ValueError(f"Expected {expected} rows in {path}; found {count}")
    # The current APDL export is in mm.  Keep the asset in metres so the
    # downstream modal projection and MATLAB path model use SI geometry.
    if np.nanmedian(np.hypot(xyz[:, 0], xyz[:, 1])) > 1.0:
        xyz *= 1e-3
        displacement *= 1e-3
    return xyz, displacement


def sensor_modes(
    path: Path,
    mode_count: int,
    nodes_per_sensor: int,
    sensor_ids: tuple[str, ...] = SENSOR_IDS,
):
    """Read labelled ANSYS patches from 9- or 10-column modal CSV rows.

    The 9-column export uses ``mode_id,sensor_id_or_angle,node_id,x,y,z,ux,uy,uz``;
    the 10-column export uses ``mode_id,frequency,angle,node_id,...``.
    Coordinates and displacements are returned in metres.
    """
    points = {name: np.empty((nodes_per_sensor, 3)) for name in sensor_ids}
    values = {
        name: np.empty((mode_count, nodes_per_sensor, 3)) for name in sensor_ids
    }
    counters = {name: np.zeros(mode_count, dtype=int) for name in sensor_ids}
    reference = _reference_frequencies(path, mode_count)
    embedded_frequency = np.full(mode_count, np.nan)
    active_mode = 0
    row_format = None
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        rows = csv.reader(f)
        next(rows, None)
        for row_index, row in enumerate(rows, start=2):
            if not row or not row[0].strip():
                continue
            if len(row) not in (9, 10):
                raise ValueError(f"Expected 9 or 10 fields in row {row_index} of {path}")
            frequency_inclusive = len(row) == 10
            if row_format is None:
                row_format = frequency_inclusive
            elif frequency_inclusive != row_format:
                raise ValueError(f"Mixed sensor row formats in {path}")
            if frequency_inclusive:
                angle = int(round(float(row[2])))
                sensor_id = f"SENSOR_{angle}"
                point_slice = slice(4, 7)
                value_slice = slice(7, 10)
            else:
                label = row[1].strip()
                if label.startswith("SENSOR_"):
                    sensor_id = label
                else:
                    angle = int(round(float(label)))
                    sensor_id = f"SENSOR_{angle}"
                point_slice = slice(3, 6)
                value_slice = slice(6, 9)
            if sensor_id not in sensor_ids:
                raise ValueError(f"Unknown sensor id {sensor_id!r}")
            mode = _positive_mode_id(row[0], path, row_index) - 1
            if mode < 0 or mode >= mode_count:
                raise ValueError(f"Mode index {mode + 1} is outside the modal summary")
            if mode != active_mode:
                if mode != active_mode + 1:
                    raise ValueError(
                        f"Sensor mode_id block sequence in {path}: expected "
                        f"{active_mode + 1} or {active_mode + 2}, found {mode + 1}"
                    )
                if any(counters[name][active_mode] != nodes_per_sensor for name in sensor_ids):
                    raise ValueError(
                        f"Incomplete sensor mode_id block {active_mode + 1} in {path}"
                    )
                active_mode = mode
            if frequency_inclusive:
                expected_frequency = (
                    reference[mode] if reference is not None else
                    embedded_frequency[mode] if math.isfinite(embedded_frequency[mode]) else None
                )
                embedded_frequency[mode] = _check_frequency(
                    row[1], expected_frequency, path, row_index
                )
            local = counters[sensor_id][mode]
            if local >= nodes_per_sensor:
                raise ValueError(f"Too many {sensor_id} nodes in mode {mode + 1}")
            point = np.asarray(row[point_slice], dtype=float)
            if mode == 0:
                points[sensor_id][local] = point
            elif not np.allclose(
                point, points[sensor_id][local], rtol=0.0, atol=1e-12
            ):
                raise ValueError(
                    f"Node ordering changed for {sensor_id} in mode {mode + 1}"
                )
            values[sensor_id][mode, local] = np.asarray(row[value_slice], dtype=float)
            counters[sensor_id][mode] += 1
    for sensor_id in sensor_ids:
        if not np.all(counters[sensor_id] == nodes_per_sensor):
            raise ValueError(
                f"{sensor_id} must contain exactly {nodes_per_sensor} nodes per mode; "
                f"found {np.unique(counters[sensor_id]).tolist()}"
            )
    # Match the ring block conversion for the current N-mm-MPa APDL export.
    all_points = np.concatenate(list(points.values()), axis=0)
    if np.nanmedian(np.hypot(all_points[:, 0], all_points[:, 1])) > 1.0:
        for sensor_id in sensor_ids:
            points[sensor_id] *= 1e-3
            values[sensor_id] *= 1e-3
    return points, values


def validate_export_qa(
    source: Path,
    *,
    expected_mode_count: int | None = None,
    expected_sensor_nodes: int | None = None,
    check_marker_freshness: bool = False,
) -> dict[str, object]:
    """Check the APDL mapping QA and completion marker when supplied.

    Historical and manually assembled exports may have neither file. A QA file
    or completion marker on its own indicates an incomplete export and fails.
    Freshness is optional because a four-sensor CSV may be augmented later.
    """
    qa_path = source / "auto_sensor_check.csv"
    marker_path = source / "UNIFIED_FE_ASSET_EXPORT_FINISHED.txt"
    if not qa_path.exists() and not marker_path.exists():
        return {"available": False}
    if not qa_path.is_file() or not marker_path.is_file():
        raise ValueError(
            f"Incomplete FE export in {source}: auto_sensor_check.csv and "
            "UNIFIED_FE_ASSET_EXPORT_FINISHED.txt must both exist"
        )

    with qa_path.open("r", encoding="utf-8", errors="ignore", newline="") as stream:
        rows = csv.reader(stream)
        header = next(rows, None)
        if header is None or [cell.strip().lower() for cell in header] != ["item", "value"]:
            raise ValueError(f"Invalid auto_sensor_check.csv header in {qa_path}")
        qa: dict[str, str] = {}
        for row_index, row in enumerate(rows, start=2):
            if not row or not any(cell.strip() for cell in row):
                continue
            if len(row) != 2 or not row[0].strip() or not row[1].strip():
                raise ValueError(f"Incomplete QA row {row_index} in {qa_path}")
            item = row[0].strip().lower()
            if item in qa:
                raise ValueError(f"Duplicate QA item {item!r} in {qa_path}")
            qa[item] = row[1].strip()

    required = (
        "ring_cx", "ring_cy", "n_s0", "n_s120", "n_s240",
        "mean120", "max120", "mean240", "max240", "unit",
    )
    missing = [item for item in required if item not in qa]
    if missing:
        raise ValueError(f"Incomplete auto_sensor_check.csv in {qa_path}: missing {missing}")
    if qa["unit"] != "mm_to_m":
        raise ValueError(f"Unexpected FE export unit {qa['unit']!r} in {qa_path}")

    numbers: dict[str, float] = {}
    for item in required[:-1]:
        try:
            numbers[item] = float(qa[item])
        except ValueError as exc:
            raise ValueError(f"Invalid QA value for {item} in {qa_path}") from exc
        if not math.isfinite(numbers[item]):
            raise ValueError(f"Nonfinite QA value for {item} in {qa_path}")
    counts = [numbers[item] for item in ("n_s0", "n_s120", "n_s240")]
    if any(not value.is_integer() or value <= 0 for value in counts) or len(set(counts)) != 1:
        raise ValueError(f"FE sensor QA node counts are invalid or unequal in {qa_path}")
    if expected_sensor_nodes is not None and int(counts[0]) != expected_sensor_nodes:
        raise ValueError(
            f"FE sensor QA has {int(counts[0])} nodes, expected {expected_sensor_nodes}"
        )
    for angle in (120, 240):
        mean = numbers[f"mean{angle}"]
        maximum = numbers[f"max{angle}"]
        if mean < 0 or maximum < mean or maximum > 2.0:
            raise ValueError(
                f"FE sensor QA mapping error for {angle} degrees exceeds "
                f"the 2 mm limit or is inconsistent in {qa_path}"
            )

    marker = marker_path.read_text(encoding="utf-8", errors="ignore")
    if not marker.lstrip().startswith("OK: unified FE modal asset exported"):
        raise ValueError(f"Invalid FE export completion marker in {marker_path}")
    matches = re.findall(r"\b(NMODES|SENSOR_0|SENSOR_120|SENSOR_240)\s*=\s*(\d+)\b", marker)
    marker_counts = dict(matches)
    if len(matches) != 4 or len(marker_counts) != 4:
        raise ValueError(f"Incomplete FE export completion marker in {marker_path}")
    if int(marker_counts["NMODES"]) <= 0:
        raise ValueError(f"Invalid mode count in {marker_path}")
    if expected_mode_count is not None and int(marker_counts["NMODES"]) != expected_mode_count:
        raise ValueError(
            f"FE export marker has {marker_counts['NMODES']} modes, "
            f"expected {expected_mode_count}"
        )
    for name, key in (("SENSOR_0", "n_s0"), ("SENSOR_120", "n_s120"),
                      ("SENSOR_240", "n_s240")):
        if int(marker_counts[name]) != int(numbers[key]):
            raise ValueError(f"FE export marker {name} disagrees with {qa_path}")

    if check_marker_freshness:
        csv_paths = [
            source / name for name in (
                "modal_summary.csv", "mode_shapes_ring_teeth.csv", "mode_shapes_sensor.csv"
            )
        ]
        if any(not path.is_file() for path in csv_paths):
            raise ValueError(f"Incomplete FE modal CSV export in {source}")
        if marker_path.stat().st_mtime + 2.0 < max(path.stat().st_mtime for path in csv_paths):
            raise ValueError(f"Stale FE export completion marker in {marker_path}")
    return {
        "available": True,
        "mode_count": int(marker_counts["NMODES"]),
        "sensor_nodes": int(counts[0]),
        "max_mapping_error_mm": max(numbers["max120"], numbers["max240"]),
    }
