"""Create thesis TVMS figures by exact two-cycle repetition of saved data.

This script is intentionally plot-only: it reads the existing one-cycle CSV
results and applies k(psi + 1) = k(psi). It neither imports nor invokes the
TVMS solver, and it performs no interpolation, filtering, or smoothing.
"""
from __future__ import annotations

import csv
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "tvms"
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".matplotlib"))

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager, ft2font
from matplotlib.font_manager import FontProperties

DPI = 600
SP_COLOR = "#1f4e79"
PR_COLOR = "#a64b00"
TIMES_FONT = FontProperties(fname=r"C:\Windows\Fonts\times.ttf", size=10)
TIMES_BOLD_FONT = FontProperties(fname=r"C:\Windows\Fonts\timesbd.ttf", size=11)
CHINESE_FONT = FontProperties(fname=r"C:\Windows\Fonts\simsun.ttc", size=10)
Y_LABEL_CN = "啮合刚度"
Y_LABEL_UNIT = "(10⁸ N/m)"


def read_total_stiffness(filename: str) -> tuple[np.ndarray, np.ndarray]:
    path = OUT / filename
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    phase = np.asarray([float(row["phase_norm"]) for row in rows])
    stiffness = np.asarray([float(row["k_total_N_per_m"]) for row in rows])
    if len(phase) < 2 or not np.isclose(phase[0], 0.0) or not np.isclose(phase[-1], 1.0):
        raise RuntimeError(f"{path} is not a complete normalized cycle [0, 1].")
    if not np.isclose(stiffness[0], stiffness[-1], rtol=1e-10, atol=0.0):
        raise RuntimeError(f"{path} does not close at the cycle endpoint.")
    return phase, stiffness


def exact_two_cycles(phase: np.ndarray, stiffness: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Concatenate saved samples exactly, retaining the endpoint only once."""
    x = np.concatenate((phase[:-1], phase + 1.0))
    y = np.concatenate((stiffness[:-1], stiffness))
    if not np.array_equal(y[: len(stiffness) - 1], y[len(stiffness) - 1 : -1]):
        raise RuntimeError("The second cycle is not an exact copy of the first.")
    return x, y


def configure_style() -> None:
    try:
        font_manager.findfont("Times New Roman", fallback_to_default=False)
    except ValueError as exc:
        raise RuntimeError("Times New Roman is required but is not installed.") from exc
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.size": 10,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def verify_explicit_unit_font() -> None:
    """Fail if the requested unit glyphs are not supplied by times.ttf."""
    font_path = Path(TIMES_FONT.get_file()).resolve()
    expected_path = Path(r"C:\Windows\Fonts\times.ttf").resolve()
    if font_path != expected_path:
        raise RuntimeError(f"Unexpected unit font file: {font_path}")
    charmap = ft2font.FT2Font(str(font_path)).get_charmap()
    missing = sorted({char for char in Y_LABEL_UNIT if not char.isspace() and ord(char) not in charmap})
    if missing:
        raise RuntimeError(f"times.ttf lacks unit-label glyphs: {missing}")
    print(f"Unit-label font file: {font_path}")
    print("Unit-label glyph coverage: PASS (including Unicode superscript 8)")


def set_mixed_ylabel(ax: plt.Axes) -> None:
    """Draw Chinese and Latin ylabel parts with explicit font files.

    Unicode superscript 8 avoids Matplotlib mathtext completely, so the full
    scale/unit expression is rasterized directly from Windows ``times.ttf``.
    """
    ax.set_ylabel("")
    common = {
        "transform": ax.transAxes,
        "rotation": 90,
        "rotation_mode": "anchor",
        "ha": "center",
        "va": "center",
        "clip_on": False,
    }
    ax.text(-0.075, 0.385, Y_LABEL_CN, fontproperties=CHINESE_FONT, **common)
    ax.text(-0.075, 0.680, Y_LABEL_UNIT, fontproperties=TIMES_FONT, **common)


def style_axis(ax: plt.Axes) -> None:
    ax.set_xlim(0.0, 2.0)
    ax.set_xticks(np.arange(0.0, 2.01, 0.25))
    ax.set_xlabel("归一化啮合相位", fontproperties=CHINESE_FONT)
    set_mixed_ylabel(ax)
    for tick in (*ax.get_xticklabels(), *ax.get_yticklabels()):
        tick.set_fontproperties(TIMES_FONT)
    ax.margins(x=0)


def save_single(filename: str, phase: np.ndarray, stiffness: np.ndarray, color: str) -> None:
    fig, ax = plt.subplots(figsize=(5.7, 3.45), dpi=DPI)
    ax.plot(phase, stiffness / 1e8, color=color, linewidth=1.35)
    style_axis(ax)
    fig.tight_layout(pad=0.55)
    fig.savefig(OUT / filename, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def save_compare(filename: str, phase: np.ndarray, sp: np.ndarray, pr: np.ndarray) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(5.9, 5.4), dpi=DPI, sharex=True)
    axes[0].plot(phase, sp / 1e8, color=SP_COLOR, linewidth=1.35)
    axes[1].plot(phase, pr / 1e8, color=PR_COLOR, linewidth=1.35)
    for label, ax in zip(("(a)", "(b)"), axes):
        ax.set_xlim(0.0, 2.0)
        ax.set_xticks(np.arange(0.0, 2.01, 0.25))
        set_mixed_ylabel(ax)
        for tick in (*ax.get_xticklabels(), *ax.get_yticklabels()):
            tick.set_fontproperties(TIMES_FONT)
        ax.margins(x=0)
        ax.text(-0.105, 1.01, label, transform=ax.transAxes,
                fontproperties=TIMES_BOLD_FONT, va="bottom", ha="left")
    axes[1].set_xlabel("归一化啮合相位", fontproperties=CHINESE_FONT)
    fig.tight_layout(pad=0.65, h_pad=0.85)
    fig.savefig(OUT / filename, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    configure_style()
    verify_explicit_unit_font()
    sp_phase, sp_one = read_total_stiffness("tvms_sp.csv")
    pr_phase, pr_one = read_total_stiffness("tvms_pr.csv")
    if not np.array_equal(sp_phase, pr_phase):
        raise RuntimeError("SP and PR phase grids are not identical.")
    phase, sp_two = exact_two_cycles(sp_phase, sp_one)
    phase_pr, pr_two = exact_two_cycles(pr_phase, pr_one)
    if not np.array_equal(phase, phase_pr):
        raise RuntimeError("SP and PR extended phase grids are not identical.")
    save_single("Fig2_8a_SP_TVMS_2cycles.png", phase, sp_two, SP_COLOR)
    save_single("Fig2_8b_PR_TVMS_2cycles.png", phase, pr_two, PR_COLOR)
    save_compare("Fig2_8_TVMS_2cycles_compare.png", phase, sp_two, pr_two)
    print(f"Exact two-cycle phase range: {phase[0]:.1f} to {phase[-1]:.1f}")
    print(f"Samples: one cycle={len(sp_one)}, two cycles={len(phase)}")
    print("Second-cycle equality: SP=PASS; PR=PASS")
    for name in ("Fig2_8a_SP_TVMS_2cycles.png", "Fig2_8b_PR_TVMS_2cycles.png", "Fig2_8_TVMS_2cycles_compare.png"):
        print(OUT / name)


if __name__ == "__main__":
    main()
