"""Create thesis-local fault-TVMS figures from saved full-recurrence CSV data.

This is a plot-only transform: samples are cyclically re-indexed around the
saved maximum stiffness-loss sample and a six-mesh-period window is selected.
Neither the original CSV data nor any TVMS/fault parameter is modified.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "results" / "tvms" / ".matplotlib_fault_local"))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties


SAMPLES_PER_CYCLE = 2000
WINDOW_CYCLES = 6
CENTER_CYCLE = 3
DPI = 600
TIMES = FontProperties(fname=r"C:\Windows\Fonts\times.ttf", size=10)
TIMES_BOLD = FontProperties(fname=r"C:\Windows\Fonts\timesbd.ttf", size=10.5)
SIMSUN = FontProperties(fname=r"C:\Windows\Fonts\simsun.ttc", size=10)


def configure_style() -> None:
    mpl.rcParams.update({
        "font.family": "Times New Roman", "font.size": 10,
        "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "legend.frameon": False,
    })


def load_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    needed = ("mesh_phase_cycles", "k_healthy_N_per_m", "k_fault_N_per_m", "delta_k_N_per_m")
    if not rows or any(key not in rows[0] for key in needed):
        raise RuntimeError(f"Unexpected fault-TVMS CSV schema: {path}")
    result = {key: np.asarray([float(row[key]) for row in rows]) for key in needed}
    if len(result["mesh_phase_cycles"]) < 2 or not np.all(np.isfinite(result["delta_k_N_per_m"])):
        raise RuntimeError(f"Invalid numeric data in {path}")
    return result


def local_window(data: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, int]:
    """Exact modular crop centered on the saved maximum stiffness loss."""
    full_count = len(data["mesh_phase_cycles"]) - 1  # terminal recurrence endpoint duplicates phase zero
    if full_count % SAMPLES_PER_CYCLE:
        raise RuntimeError("Input grid is not an integer number of 2000-sample mesh periods.")
    center_index = int(np.argmax(data["delta_k_N_per_m"][:full_count]))
    if data["delta_k_N_per_m"][center_index] <= 0:
        raise RuntimeError("No fault stiffness-loss event is present in the input CSV.")
    local_indices = (np.arange(WINDOW_CYCLES * SAMPLES_PER_CYCLE + 1)
                     + center_index - CENTER_CYCLE * SAMPLES_PER_CYCLE) % full_count
    x = np.arange(WINDOW_CYCLES * SAMPLES_PER_CYCLE + 1) / SAMPLES_PER_CYCLE
    healthy_curve = data["k_healthy_N_per_m"][local_indices]
    fault_curve = data["k_fault_N_per_m"][local_indices]
    if not np.isclose(x[CENTER_CYCLE * SAMPLES_PER_CYCLE], CENTER_CYCLE):
        raise RuntimeError("Local figure center is not exactly N=3.")
    if not np.isclose(healthy_curve[CENTER_CYCLE * SAMPLES_PER_CYCLE], data["k_healthy_N_per_m"][center_index]):
        raise RuntimeError("Cyclic extraction changed a healthy sample.")
    if not np.isclose(fault_curve[CENTER_CYCLE * SAMPLES_PER_CYCLE], data["k_fault_N_per_m"][center_index]):
        raise RuntimeError("Cyclic extraction changed a fault sample.")
    return x, healthy_curve, fault_curve, center_index / SAMPLES_PER_CYCLE, full_count // SAMPLES_PER_CYCLE


def set_mixed_ylabel(ax: plt.Axes, chinese: str) -> None:
    """Keep Chinese in SimSun and full unit expression directly in times.ttf."""
    ax.set_ylabel("")
    common = {"transform": ax.transAxes, "rotation": 90, "rotation_mode": "anchor",
              "ha": "center", "va": "center", "clip_on": False}
    ax.text(-0.095, 0.36, chinese, fontproperties=SIMSUN, **common)
    ax.text(-0.095, 0.73, "(10⁸ N/m)", fontproperties=TIMES, **common)


def style_axis(ax: plt.Axes, ylabel: str, panel: str | None = None) -> None:
    ax.set_xlim(0, WINDOW_CYCLES)
    ax.set_xticks(np.arange(0, WINDOW_CYCLES + 1, 1))
    set_mixed_ylabel(ax, ylabel)
    ax.axvline(CENTER_CYCLE, color="#8e8e8e", linestyle="--", linewidth=0.65, zorder=0)
    if panel:
        ax.text(-0.10, 1.01, panel, transform=ax.transAxes, fontproperties=TIMES_BOLD,
                ha="left", va="bottom")
    for text in (*ax.get_xticklabels(), *ax.get_yticklabels()):
        text.set_fontproperties(TIMES)
    ax.margins(x=0)


def plot_single(path: Path, data: dict[str, np.ndarray], ylabel: str) -> tuple[float, int]:
    x, healthy_curve, fault_curve, source_center, recurrence = local_window(data)
    fig, ax = plt.subplots(figsize=(5.8, 3.45), dpi=DPI)
    ax.plot(x, healthy_curve / 1e8, color="#000000", linewidth=0.80, label="健康")
    ax.plot(x, fault_curve / 1e8, color="#1f4e79", linewidth=1.35, label="故障")
    style_axis(ax, ylabel)
    ax.set_xlabel("啮合相位（周期）", fontproperties=SIMSUN)
    ax.legend(prop=SIMSUN, loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8.5)
    fig.subplots_adjust(left=0.18, right=0.78, bottom=0.18, top=0.95)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return source_center, recurrence


def plot_planet(path: Path, sp_data: dict[str, np.ndarray], pr_data: dict[str, np.ndarray]) -> tuple[float, float, int]:
    sp_x, sp_h, sp_f, sp_center, sp_recurrence = local_window(sp_data)
    pr_x, pr_h, pr_f, pr_center, pr_recurrence = local_window(pr_data)
    if sp_recurrence != pr_recurrence:
        raise RuntimeError("SP and PR planet-fault recurrence lengths differ.")
    fig, axes = plt.subplots(2, 1, figsize=(5.95, 5.25), dpi=DPI, sharex=True)
    for ax, x, h, f, label, color, panel in (
        (axes[0], sp_x, sp_h, sp_f, "SP 啮合刚度", "#1f4e79", "(a)"),
        (axes[1], pr_x, pr_h, pr_f, "PR 啮合刚度", "#a64b00", "(b)"),
    ):
        ax.plot(x, h / 1e8, color="#000000", linewidth=0.80, label="健康")
        ax.plot(x, f / 1e8, color=color, linewidth=1.35, label="故障")
        style_axis(ax, label, panel)
        ax.legend(prop=SIMSUN, loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8.5)
    axes[1].set_xlabel("啮合相位（周期）", fontproperties=SIMSUN)
    fig.subplots_adjust(left=0.18, right=0.78, bottom=0.13, top=0.94, hspace=0.16)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return sp_center, pr_center, sp_recurrence


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path, help="Existing full-recurrence local-fault TVMS results directory")
    parser.add_argument("--q-tag", default="q050", choices=("q025", "q050"),
                        help="Existing q case used for the two formal paper figures (default: q050).")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    sun_source = run_dir / f"sun_fault_sp_p1_{args.q_tag}.csv"
    planet_sp_source = run_dir / f"planet_fault_sp_p1_{args.q_tag}.csv"
    planet_pr_source = run_dir / f"planet_fault_pr_p1_{args.q_tag}.csv"
    for source in (sun_source, planet_sp_source, planet_pr_source):
        if not source.is_file():
            raise FileNotFoundError(source)
    configure_style()
    sun_center, sun_recurrence = plot_single(run_dir / "Fig2_9_sun_fault_TVMS_local.png",
                                              load_csv(sun_source), "SP 啮合刚度")
    sp_center, pr_center, planet_recurrence = plot_planet(
        run_dir / "Fig2_10_planet_fault_TVMS_local.png",
        load_csv(planet_sp_source), load_csv(planet_pr_source))
    physical_forward_interval = (pr_center - sp_center) % planet_recurrence
    if not 15.0 < physical_forward_interval < 16.0:
        raise RuntimeError(f"Unexpected preserved SP->PR interval: {physical_forward_interval:.6f} cycles")
    summary = [
        "Local fault-TVMS paper-figure transform", "",
        "Scope: plot-only. Full-recurrence fault CSV files were read without modification.",
        f"Selected existing severity case: {args.q_tag}.",
        "Transform: exact cyclic sample re-indexing and a [0,6]-cycle local crop; no interpolation, smoothing or TVMS recomputation.",
        "At N=3 the saved maximum local stiffness-loss sample is placed at the figure center.",
        f"Sun SP source center={sun_center:.6f} displayed cycles; raw recurrence={sun_recurrence} cycles.",
        f"Planet SP source center={sp_center:.6f} displayed cycles; planet PR source center={pr_center:.6f} displayed cycles; raw recurrence={planet_recurrence} cycles.",
        f"Preserved forward physical SP->PR event interval: delta_N={physical_forward_interval:.6f} mesh cycles (approximately 15.5 cycles).",
        "SP and PR are separately re-centered only for the local comparison figure; their original full-range CSV time relation is unchanged.",
        "",
        "Input CSV SHA-256:",
        f"{sun_source.name}: {sha256(sun_source)}",
        f"{planet_sp_source.name}: {sha256(planet_sp_source)}",
        f"{planet_pr_source.name}: {sha256(planet_pr_source)}",
    ]
    (run_dir / "local_fault_tvms_plot_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))
    print(run_dir / "Fig2_9_sun_fault_TVMS_local.png")
    print(run_dir / "Fig2_10_planet_fault_TVMS_local.png")


if __name__ == "__main__":
    main()
