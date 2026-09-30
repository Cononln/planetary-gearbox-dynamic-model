#!/usr/bin/env python3
"""Reconstruct model-relative ring-tooth labels for the 48 saved fault events.

The labels are a bookkeeping convention, not surveyed ANSYS/CAD tooth IDs.
R84/R1 straddles zero of the source model's local carrier/planet angle.  At a
phase between tooth pitches, the two bracketing labels describe possible PR
contact teeth; they do not claim that both teeth carry an equal load.

Inputs are the saved corrected event table and the *same* fault-source MAT
whose carrierAngleRad was used for the FE re-registration.  Event times are
relative to that retained MAT record (after its 0.2 s discarded transient).
No dynamics or FE response is recomputed.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "fe_causal_q090_corrected_registration_20260928"
DEFAULT_EVENTS = RESULTS / "corrected_event_metrics.csv"
DEFAULT_METADATA = RESULTS / "README.json"
DEFAULT_OUTPUT = RESULTS / "mapping_evidence" / "model_relative_ring_tooth_sequence.csv"
RING_TEETH = 84
PLANET_PHASE_RAD = (0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0)
EXPECTED_FIRST_FOUR_PAIRS = ("R63/R64", "R42/R43", "R21/R22", "R84/R1")


def label(tooth_index: int) -> str:
    """Map the modulo-84 index 0 to R84 and 1..83 to R1..R83."""

    return f"R{tooth_index % RING_TEETH or RING_TEETH}"


def ring_position(carrier_angle_rad: float, planet: int) -> float:
    """Unwrapped carrier-plus-planet angle measured in ring-tooth pitches."""

    if planet not in (1, 2, 3):
        raise ValueError(f"Invalid active_planet={planet}; expected 1, 2, or 3")
    return RING_TEETH * (carrier_angle_rad + PLANET_PHASE_RAD[planet - 1]) / (2.0 * math.pi)


def bracketing_teeth(position: float) -> tuple[str, str]:
    lower = math.floor(position)
    return label(lower), label(lower + 1)


def read_events(path: Path) -> list[dict[str, str]]:
    """Collapse four sensor rows per saved event, checking source fields agree."""

    source_keys = (
        "active_planet",
        "source_event_start_time_s",
        "source_event_peak_time_s",
        "source_event_end_time_s",
        "contact_angle_global_deg",
    )
    grouped: dict[int, list[dict[str, str]]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["registration"] != "corrected":
                raise ValueError(f"Expected corrected registration in {path}")
            grouped.setdefault(int(row["event_number"]), []).append(row)
    if sorted(grouped) != list(range(1, 49)):
        raise ValueError(f"Expected saved events 1..48; found {sorted(grouped)}")
    events: list[dict[str, str]] = []
    for event_number in range(1, 49):
        group = grouped[event_number]
        if len(group) != 4 or {float(row["sensor_angle_local_deg"]) for row in group} != {0.0, 90.0, 120.0, 240.0}:
            raise ValueError(f"Event {event_number} lacks the four saved sensor rows")
        first = group[0]
        if any(any(row[key] != first[key] for key in source_keys) for row in group[1:]):
            raise ValueError(f"Event {event_number} has inconsistent source fields")
        sensor_zero = next(row for row in group if float(row["sensor_angle_local_deg"]) == 0.0)
        events.append({**first, "sensor0_global_deg": sensor_zero["sensor_angle_global_deg"]})
    return events


def source_samples(path: Path, n_samples: int) -> tuple[np.ndarray, np.ndarray]:
    """Read retained time and carrier-angle columns from a MATLAB v7.3 MAT."""

    try:
        import h5py
    except ImportError as exc:
        raise RuntimeError("Reading the source MATLAB v7.3 MAT requires h5py") from exc
    with h5py.File(path, "r") as handle:
        for key in ("time", "carrierAngleRad"):
            if key not in handle or handle[key].shape != (1, handle[key].shape[-1]):
                raise ValueError(f"Expected MATLAB row vector {key} in {path}")
            if handle[key].shape[-1] < n_samples:
                raise ValueError(f"{key} has fewer than {n_samples} saved samples")
        time = np.asarray(handle["time"][0, :n_samples], dtype=float)
        carrier = np.asarray(handle["carrierAngleRad"][0, :n_samples], dtype=float)
    if not (np.all(np.isfinite(time)) and np.all(np.isfinite(carrier)) and np.all(np.diff(time) > 0)):
        raise ValueError("Source time/carrierAngleRad must be finite and time increasing")
    return time, carrier


def sample_index(time: np.ndarray, event_time: float) -> int:
    """Match an event-table clock value to the exact retained MAT sample."""

    right = int(np.searchsorted(time, event_time))
    candidates = [i for i in (right - 1, right) if 0 <= i < time.size]
    index = min(candidates, key=lambda i: abs(float(time[i]) - event_time))
    if abs(float(time[index]) - event_time) > 1e-10:
        raise ValueError(f"Event time {event_time:.12g} is not a retained MAT sample")
    return index


def reconstruct(events: list[dict[str, str]], time: np.ndarray, carrier: np.ndarray) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for event_number, event in enumerate(events, start=1):
        planet = int(event["active_planet"])
        start_time = float(event["source_event_start_time_s"])
        peak_time = float(event["source_event_peak_time_s"])
        end_time = float(event["source_event_end_time_s"])
        start, peak, end = (sample_index(time, t) for t in (start_time, peak_time, end_time))
        if not start <= peak <= end:
            raise ValueError(f"Event {event_number} has invalid start/peak/end order")

        positions = [ring_position(float(angle), planet) for angle in carrier[start : end + 1]]
        if any(b <= a for a, b in zip(positions, positions[1:])):
            raise ValueError(f"Event {event_number} carrier phase is not increasing")
        peak_position = ring_position(float(carrier[peak]), planet)
        peak_pair = bracketing_teeth(peak_position)
        # Consecutive samples traverse the full interval here.  Include both
        # bracketing labels at every retained sample in source-time order.
        whole_event_teeth = list(dict.fromkeys(
            tooth for position in positions for tooth in bracketing_teeth(position)
        ))
        if len(whole_event_teeth) != 3:
            raise ValueError(
                f"Event {event_number} spans {len(whole_event_teeth)} labels, "
                f"expected three: {whole_event_teeth}"
            )

        # The corrected FE contact angle is reflected around SENSOR_0's
        # exported global direction.  It is an independent check that the
        # MAT carrier angle is the one used to create this event table.
        predicted_contact = (float(event["sensor0_global_deg"]) - math.degrees(
            float(carrier[peak]) + PLANET_PHASE_RAD[planet - 1]
        )) % 360.0
        saved_contact = float(event["contact_angle_global_deg"])
        angular_error = (predicted_contact - saved_contact + 180.0) % 360.0 - 180.0
        if abs(angular_error) > 1e-6:
            raise ValueError(f"Event {event_number} MAT/FE contact-angle mismatch: {angular_error:.3g} deg")

        output.append({
            "event_number": event_number,
            "active_planet": planet,
            "source_event_start_time_s": start_time,
            "source_event_peak_time_s": peak_time,
            "source_event_end_time_s": end_time,
            "carrier_angle_peak_rad": float(carrier[peak]),
            "peak_model_ring_phase_tooth_pitches": peak_position,
            "peak_ring_tooth_pair": "/".join(peak_pair),
            "whole_event_ring_tooth_set": ";".join(whole_event_teeth),
            "whole_event_tooth_count": len(whole_event_teeth),
        })
    actual_first_four = tuple(str(row["peak_ring_tooth_pair"]) for row in output[:4])
    if actual_first_four != EXPECTED_FIRST_FOUR_PAIRS:
        raise ValueError(f"First four peak pairs {actual_first_four} differ from {EXPECTED_FIRST_FOUR_PAIRS}")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events", type=Path, default=DEFAULT_EVENTS)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--source-mat", type=Path, help="Override source_fault in README.json")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    metadata = json.loads(args.metadata.read_text(encoding="utf-8"))
    source_path = args.source_mat or Path(metadata["source_fault"])
    events = read_events(args.events)
    time, carrier = source_samples(source_path, int(metadata["n_samples"]))
    rows = reconstruct(events, time, carrier)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} model-relative ring-tooth events to {args.out}")


if __name__ == "__main__":
    main()
