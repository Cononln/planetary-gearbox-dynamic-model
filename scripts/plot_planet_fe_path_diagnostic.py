#!/usr/bin/env python3
"""Plot one exploratory figure from the auditable planet-crack event CSV.

Claim: the strongest radial sensor at a planet-crack event depends on FE
modal projection as well as PR contact angle; nearest is not always strongest.
Every point is one complete SP or PR event from the same q_c=0.90 mm source.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "results" / "fe_planet_q090_23mode_exploratory_20260929"
sys.path.insert(0, str(PROJECT / "scripts"))
from figure_qa import require_matplotlib_panel_alignment  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reflected-normal", action="store_true", help="Plot the geometrically reflected PR-normal replay without replacing the original diagnostic.")
    args = parser.parse_args()
    data_dir = SOURCE / "reflected_normal_sensitivity" if args.reflected_normal else SOURCE
    rows = list(csv.DictReader((data_dir / "planet_sensor_event_metrics.csv").open(encoding="utf-8-sig", newline="")))
    if len(rows) != 40:
        raise ValueError(f"Expected ten complete events at four sensors, got {len(rows)} rows.")
    mpl.rcParams.update({
        "font.family": ["Times New Roman", "SimSun", "Arial"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 8.5,
        "axes.labelsize": 8.5,
        "axes.titlesize": 9,
        "axes.linewidth": 0.65,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })
    cn = FontProperties(family="SimSun", size=8.5)
    cn_title = FontProperties(family="SimSun", size=9)
    colors = {0: "#1f77b4", 90: "#c44e52", 120: "#228c67", 240: "#7a5aa6"}
    shapes = {"PR": "o", "SP": "^"}
    fig, ax = plt.subplots(figsize=(7.0, 4.1), constrained_layout=False)
    for sensor, color in colors.items():
        for kind, marker in shapes.items():
            data = [r for r in rows if int(float(r["sensor_local_deg"])) == sensor and r["source_mesh"] == kind]
            ax.scatter(
                [float(r["pr_contact_angle_local_deg"]) for r in data],
                [float(r["peak_per_pr1_force_delta_ms2_per_N"]) * 1000.0 for r in data],
                marker=marker, s=43, color=color, edgecolor="black", linewidth=0.35,
                zorder=3,
            )
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    peak_max = max(float(row["peak_per_pr1_force_delta_ms2_per_N"]) * 1000.0 for row in rows)
    ax.set_ylim(0, peak_max * 1.12 if args.reflected_normal else 5.35)
    ax.set_xlabel("PR 接触位置 / °（自 0° 测点顺时针）", fontproperties=cn, labelpad=7)
    ax.set_ylabel("Event peak ratio / [10⁻³ (m/s²)/N]", labelpad=7)
    ax.set_title("行星轮裂纹：故障事件位置与四测点响应", fontproperties=cn_title, pad=39)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    sensor_handles = [Line2D([], [], marker="o", linestyle="None", color=color, markersize=5.5, label=f"Sensor {sensor}°") for sensor, color in colors.items()]
    mesh_handles = [Line2D([], [], marker=shapes[kind], linestyle="None", color="#333333", markersize=5.5, label=f"{kind} event") for kind in ("SP", "PR")]
    ax.legend(handles=sensor_handles + mesh_handles, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.01), frameon=False, handletextpad=0.35, columnspacing=1.2)
    note = "23 阶可采样模态；5 次 SP + 5 次 PR 事件；镜像法向；纵轴为故障增量峰值/P1-PR力增量峰值。" if args.reflected_normal else "23 阶可采样模态；5 次 SP + 5 次 PR 完整事件；幅值为故障减健康后的峰值比。"
    fig.text(0.12, 0.025, note, fontproperties=FontProperties(family="SimSun", size=8), color="#444444")
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.20, top=0.77)

    stem = "Fig_planet_q090_FE_path_event_angle_reflected_exploratory" if args.reflected_normal else "Fig_planet_q090_FE_path_event_angle_exploratory"
    base = SOURCE / stem
    require_matplotlib_panel_alignment(
        fig,
        json_out=str(base) + ".alignment.json",
        overlay_svg=str(base) + ".alignment.svg",
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        strict=True,
    )
    fig.savefig(str(base) + ".pdf", bbox_inches="tight")
    fig.savefig(str(base) + ".svg", bbox_inches="tight")
    fig.savefig(str(base) + ".png", dpi=600, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
