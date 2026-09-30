#!/usr/bin/env python3
"""Compare the two separately computed PR-normal hypotheses by event."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1] / "results" / "fe_planet_q090_23mode_exploratory_20260929"
ALTERNATE = ROOT / "reflected_normal_sensitivity"


def read(path: Path) -> dict[int, dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return {int(row["event_number"]): row for row in csv.DictReader(handle)}


def main() -> None:
    original_events = read(ROOT / "planet_fault_events.csv")
    alternate_events = read(ALTERNATE / "planet_fault_events.csv")
    if original_events.keys() != alternate_events.keys():
        raise ValueError("The two FE hypotheses have different source events.")
    with (ROOT / "planet_normal_hypothesis_comparison.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        fields = ["event_number", "source_mesh", "contact_local_deg", "nearest_sensor_deg", "original_strongest_sensor_deg", "reflected_strongest_sensor_deg", "same_strongest_sensor"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for number in original_events:
            a, b = original_events[number], alternate_events[number]
            for key in ("source_mesh", "pr_contact_angle_local_deg_at_source_peak", "nearest_sensor_local_deg"):
                if a[key] != b[key]:
                    raise ValueError(f"Event {number} differs in {key}.")
            writer.writerow({
                "event_number": number,
                "source_mesh": a["source_mesh"],
                "contact_local_deg": a["pr_contact_angle_local_deg_at_source_peak"],
                "nearest_sensor_deg": a["nearest_sensor_local_deg"],
                "original_strongest_sensor_deg": a["strongest_sensor_local_deg"],
                "reflected_strongest_sensor_deg": b["strongest_sensor_local_deg"],
                "same_strongest_sensor": int(a["strongest_sensor_local_deg"] == b["strongest_sensor_local_deg"]),
            })
    def amplitudes(path: Path) -> dict[tuple[int, float], float]:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            return {(int(row["event_number"]), float(row["sensor_local_deg"])): float(row["fault_increment_peak_abs_ms2"]) for row in csv.DictReader(handle)}
    orig = amplitudes(ROOT / "planet_sensor_event_metrics.csv")
    alt = amplitudes(ALTERNATE / "planet_sensor_event_metrics.csv")
    ratios = np.array([alt[key] / orig[key] for key in orig])
    print(f"same four-sensor event winner: {sum(original_events[i]['strongest_sensor_local_deg'] == alternate_events[i]['strongest_sensor_local_deg'] for i in original_events)}/{len(original_events)}")
    print(f"amplitude ratio alternate/original across 40 event-sensor pairs: min={ratios.min():.4f}, median={np.median(ratios):.4f}, max={ratios.max():.4f}")


if __name__ == "__main__":
    main()
