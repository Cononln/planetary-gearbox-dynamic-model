"""Plot-only thesis figures for Section 2.3.2 root-crack TVMS results.

The script reads saved healthy/fault CSV data and applies the established exact
cyclic local-window transform.  It does not evaluate or modify the TVMS model.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "results" / "tvms" / ".matplotlib_section_232"))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties

from plot_fault_tvms_local import CENTER_CYCLE, DPI, configure_style, load_csv, local_window


OUT = ROOT / "results" / "tvms"
CM_PER_INCH = 2.54
TIMES = FontProperties(fname=r"C:\Windows\Fonts\times.ttf", size=11)
TIMES_BOLD = FontProperties(fname=r"C:\Windows\Fonts\timesbd.ttf", size=11)
SIMSUN = FontProperties(fname=r"C:\Windows\Fonts\simsun.ttc", size=11)
SIMSUN_SMALL = FontProperties(fname=r"C:\Windows\Fonts\simsun.ttc", size=10)
HEALTHY = "#000000"
SUN_FAULT = "#b1292f"
PLANET_FAULT = "#1f5f91"


def setup() -> None:
    configure_style()
    mpl.rcParams.update({
        "font.family": "Times New Roman", "font.size": 11,
        "mathtext.fontset": "custom", "mathtext.rm": "Times New Roman",
        "mathtext.it": "Times New Roman:italic", "mathtext.bf": "Times New Roman:bold",
        "axes.linewidth": 0.8, "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "svg.fonttype": "none", "pdf.fonttype": 42,
        "axes.formatter.useoffset": False,
    })


def add_mixed_ylabel(fig: plt.Figure, ax: plt.Axes, mesh: str) -> None:
    """Single compact label; all ticks/variables retain Times New Roman."""
    ax.set_ylabel(rf"{mesh} 啮合刚度（$10^8\ \mathrm{{N/m}}$）", fontproperties=SIMSUN, labelpad=8)


def style_axis(ax: plt.Axes, fig: plt.Figure, mesh: str, panel: str | None = None) -> None:
    ax.set_xlim(0, 6)
    ax.set_xticks(np.arange(0, 7, 1))
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax.yaxis.get_offset_text().set_visible(False)
    ax.axvline(CENTER_CYCLE, color="#a6a6a6", linestyle="--", linewidth=0.65, zorder=0)
    ax.tick_params(axis="both", direction="out", width=0.75, length=3.8, pad=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for tick in (*ax.get_xticklabels(), *ax.get_yticklabels()):
        tick.set_fontproperties(TIMES)
    if panel:
        ax.text(0.01, 0.97, panel, transform=ax.transAxes, fontproperties=TIMES_BOLD,
                va="top", ha="left")
    add_mixed_ylabel(fig, ax, mesh)


def draw(ax: plt.Axes, x: np.ndarray, health: np.ndarray, fault: np.ndarray, color: str) -> tuple[object, object]:
    healthy_line, = ax.plot(x, health / 1e8, color=HEALTHY, linewidth=0.85, label="健康")
    fault_line, = ax.plot(x, fault / 1e8, color=color, linewidth=1.35, label="裂纹故障")
    return healthy_line, fault_line


def save(fig: plt.Figure, png: Path) -> None:
    fig.savefig(png, dpi=DPI)
    fig.savefig(png.with_suffix(".svg"))
    fig.savefig(png.with_suffix(".pdf"))
    plt.close(fig)


def plot_sun(source: Path) -> None:
    x, health, fault, _, _ = local_window(load_csv(source))
    fig, ax = plt.subplots(figsize=(14 / CM_PER_INCH, 8 / CM_PER_INCH), dpi=DPI)
    handles = draw(ax, x, health, fault, SUN_FAULT)
    style_axis(ax, fig, "SP")
    ax.set_xlabel("啮合相位（周期）", fontproperties=SIMSUN, labelpad=5)
    fig.legend(handles, ("健康", "裂纹故障"), prop=SIMSUN_SMALL, frameon=False, ncol=2,
               loc="upper center", bbox_to_anchor=(0.55, 0.985), handlelength=2.1,
               handletextpad=0.5, columnspacing=1.4)
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.20, top=0.85)
    save(fig, OUT / "Fig2_10_sun_root_crack_TVMS.png")


def plot_planet(sp_source: Path, pr_source: Path) -> None:
    sp_x, sp_h, sp_f, _, _ = local_window(load_csv(sp_source))
    pr_x, pr_h, pr_f, _, _ = local_window(load_csv(pr_source))
    fig, axes = plt.subplots(2, 1, figsize=(14 / CM_PER_INCH, 13 / CM_PER_INCH), dpi=DPI, sharex=True)
    handles = draw(axes[0], sp_x, sp_h, sp_f, PLANET_FAULT)
    draw(axes[1], pr_x, pr_h, pr_f, PLANET_FAULT)
    style_axis(axes[0], fig, "SP", "(a)")
    style_axis(axes[1], fig, "PR", "(b)")
    axes[1].set_xlabel("啮合相位（周期）", fontproperties=SIMSUN, labelpad=5)
    fig.legend(handles, ("健康", "裂纹故障"), prop=SIMSUN_SMALL, frameon=False, ncol=2,
               loc="upper center", bbox_to_anchor=(0.55, 0.992), handlelength=2.1,
               handletextpad=0.5, columnspacing=1.4)
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.12, top=0.90, hspace=0.18)
    save(fig, OUT / "Fig2_11_planet_root_crack_TVMS.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path,
                        default=OUT / "crack_tvms_pengyue_qc0p900mm_rerun",
                        help="Saved full-recurrence root-crack-TVMS result directory.")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    names = ("sun_fault_sp_p1_qc0p900mm.csv", "planet_fault_sp_p1_qc0p900mm.csv", "planet_fault_pr_p1_qc0p900mm.csv")
    sources = [run_dir / name for name in names]
    missing = [str(source) for source in sources if not source.is_file()]
    if missing:
        raise FileNotFoundError("Missing source CSV:\n" + "\n".join(missing))
    OUT.mkdir(parents=True, exist_ok=True)
    setup()
    plot_sun(sources[0])
    plot_planet(sources[1], sources[2])
    print("Plot-only exports complete: Fig2_10 and Fig2_11.")


if __name__ == "__main__":
    main()
