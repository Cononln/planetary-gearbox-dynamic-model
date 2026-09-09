"""Replot saved local-fault TVMS CSV results; no stiffness calculation."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from run_fault_tvms import plot_planet, plot_sun, configure_plot


def read_rows(path: Path) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            row: dict[str, float | int] = {}
            for key, value in raw.items():
                if value == "":
                    row[key] = ""
                elif key in {"mesh_cycle_index", "active_pair_count", "fault_active_pair_count",
                             "entry_pair_1", "entry_pair_2", "first_tooth_pair_1", "first_tooth_pair_2",
                             "second_tooth_pair_1", "second_tooth_pair_2", "is_fault_pair_1", "is_fault_pair_2"}:
                    row[key] = int(value)
                else:
                    row[key] = float(value)
                # plot functions only need numeric stiffness/phase fields;
                # retaining all fields preserves the original row schema.
            rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path, help="Existing fault_tvms result directory")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    configure_plot()
    for tag in ("q025", "q050"):
        q = int(tag[1:]) / 100
        sun = {planet: read_rows(run_dir / f"sun_fault_sp_p{planet}_{tag}.csv") for planet in (1, 2, 3)}
        sp = read_rows(run_dir / f"planet_fault_sp_p1_{tag}.csv")
        pr = read_rows(run_dir / f"planet_fault_pr_p1_{tag}.csv")
        plot_sun(run_dir / f"Fig_fault_sun_{tag}.png", q, sun)
        plot_planet(run_dir / f"Fig_fault_planet_{tag}.png", q, sp, pr)
    print(f"Replotted saved CSV data only: {run_dir}")


if __name__ == "__main__":
    main()
