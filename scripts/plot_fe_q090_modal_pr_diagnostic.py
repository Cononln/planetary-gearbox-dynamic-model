#!/usr/bin/env python3
"""Plot data-backed diagnostics for the 0.90 mm sun-crack FE path.

The plots use the saved modal/PR-input decomposition. They do not alter the
18-DOF source response, the ANSYS modes, or any sensor trace.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = DEFAULT_PROJECT / "results" / "fe_q090_modal_pr_decomposition_20260929"
sys.path.insert(0, str(DEFAULT_PROJECT / "scripts"))
BLUE = "#245B9A"
ORANGE = "#CF6A35"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def configure() -> None:
    mpl.rcParams.update({
        "font.family": ["Times New Roman", "SimSun", "sans-serif"],
        "font.size": 10,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    })


def save_diagnostic(fig: plt.Figure, stem: Path) -> None:
    fig.canvas.draw()
    from figure_qa import require_matplotlib_panel_alignment

    require_matplotlib_panel_alignment(
        fig,
        json_out=str(stem) + ".alignment.json",
        overlay_svg=str(stem) + ".alignment.svg",
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        strict=True,
    )
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def plot_four_states(summary: list[dict[str, str]], out_dir: Path) -> None:
    by_key = {(int(float(row["state_centre_local_deg"])), int(row["sensor_local_deg"])): row
              for row in summary}
    angles = [0, 90, 180, 270]
    s0 = np.array([float(by_key[a, 0]["mean_total_delta_abs_peak_ms2"]) for a in angles]) * 1e4
    s90 = np.array([float(by_key[a, 90]["mean_total_delta_abs_peak_ms2"]) for a in angles]) * 1e4
    if len(summary) != 8 or any(int(row["event_count"]) != 12 for row in summary):
        raise ValueError("Expected four contact states, two sensors, and 12 events per state.")

    fig, ax = plt.subplots(figsize=(6.7, 3.7), layout="constrained")
    x = np.arange(4)
    width = 0.34
    ax.bar(x - width / 2, s0, width, color=BLUE, label="0°测点")
    ax.bar(x + width / 2, s90, width, color=ORANGE, label="90°测点")
    ax.set_xticks(x, [f"约{a}°" for a in angles])
    ax.set_xlabel("PR 接触方位（相对0°测点，顺时针）")
    ax.set_ylabel("故障增量峰值（10⁻⁴ m/s²）")
    ax.set_title("太阳轮裂纹：四类接触方位的测点响应", pad=10)
    ax.legend(ncol=2, loc="upper right")
    ax.set_ylim(0, max(np.max(s0), np.max(s90)) * 1.16)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.45)
    ax.set_axisbelow(True)
    save_diagnostic(fig, out_dir / "sensor_peak_by_pr_contact")


def plot_counterexample_modes(rows: list[dict[str, str]], out_dir: Path) -> None:
    # State 1 is the counterexample: local PR contact ~272.8 degrees. It is
    # angularly closer to SENSOR_0, while SENSOR_90 has the larger peak.
    selected = [row for row in rows if int(row["state"]) == 1]
    if len(selected) != 12 * 2 * 60:
        raise ValueError("Expected all 60 modes for 12 events and two sensors in state 1.")
    shares: dict[int, dict[int, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in selected:
        sensor = int(row["sensor_local_deg"])
        mode = int(row["mode"])
        shares[sensor][mode].append(float(row["projection_share_event_window"]))
    averages = {sensor: {mode: float(np.mean(values)) for mode, values in modes.items()}
                for sensor, modes in shares.items()}
    if set(averages) != {0, 90} or any(len(averages[sensor]) != 60 for sensor in (0, 90)):
        raise ValueError("Modal contribution table is incomplete.")
    top_s0 = sorted(averages[0], key=lambda m: -abs(averages[0][m]))[:5]
    top_s90 = sorted(averages[90], key=lambda m: -abs(averages[90][m]))[:5]
    modes = sorted(set(top_s0 + top_s90), key=lambda m: -max(abs(averages[0][m]), abs(averages[90][m])))
    labels = [f"M{mode}" for mode in modes] + ["其余模态"]
    values = {}
    for sensor in (0, 90):
        shown = [averages[sensor][mode] for mode in modes]
        values[sensor] = np.array(shown + [1.0 - sum(shown)])
        if not np.isclose(sum(averages[sensor].values()), 1.0, atol=1e-8):
            raise ValueError(f"Signed modal projections fail to sum to one for sensor {sensor}.")

    fig, ax = plt.subplots(figsize=(7.3, 4.5), layout="constrained")
    x = np.arange(len(labels))
    width = 0.35
    ax.bar(x - width / 2, values[0], width, color=BLUE, label="0°测点")
    ax.bar(x + width / 2, values[90], width, color=ORANGE, label="90°测点")
    ax.axhline(0, color="#333333", linewidth=0.75)
    ax.set_xticks(x, labels)
    ax.set_ylabel("事件窗有符号模态投影份额")
    ax.set_xlabel("模态编号")
    ax.set_title("PR接触约272.8°：两测点的模态贡献", pad=10)
    ax.legend(ncol=2, loc="upper right")
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.45)
    ax.set_axisbelow(True)
    save_diagnostic(fig, out_dir / "counterexample_signed_modal_projection")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_DATA)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    configure()
    plot_four_states(read_rows(args.data_dir / "state_summary.csv"), args.out_dir)
    plot_counterexample_modes(read_rows(args.data_dir / "event_mode_sums.csv"), args.out_dir)
    print(f"Saved diagnostic figures to {args.out_dir}")


if __name__ == "__main__":
    main()
