"""Export separate 14 cm x 8 cm fault-TVMS thesis figures from saved CSVs."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "results" / "tvms" / ".matplotlib_fault_compact"))

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.font_manager import FontProperties

from plot_fault_tvms_local import CENTER_CYCLE, DPI, configure_style, load_csv, local_window


CM_PER_INCH = 2.54
TIMES = FontProperties(fname=r"C:\Windows\Fonts\times.ttf", size=11)
SIMSUN = FontProperties(fname=r"C:\Windows\Fonts\simsun.ttc", size=11)


def add_mixed_ylabel(fig: plt.Figure, mesh_label: str) -> None:
    """Place a single vertical label with SimSun Chinese and Times unit text."""
    common = {
        "transform": fig.transFigure,
        "rotation": 90,
        "rotation_mode": "anchor",
        "ha": "center",
        "va": "center",
        "clip_on": False,
    }
    # Each language run is an independent Text artist.  Matplotlib otherwise
    # applies the Chinese fallback font to a mixed-language label, including
    # the PR/SP identifier and unit string.
    parts = (
        (mesh_label, TIMES),
        (" 啮合刚度（", SIMSUN),
        ("10⁸ N/m", TIMES),
        ("）", SIMSUN),
    )
    artists = [fig.text(0.075, 0.5, text, fontproperties=font, **common)
               for text, font in parts]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    heights = [artist.get_window_extent(renderer).height for artist in artists]
    gap = fig.dpi * 1.0 / 72.0  # one point: clear separation without visual breaks
    total_height = sum(heights) + gap * (len(artists) - 1)
    figure_height = fig.bbox.height
    lower = (figure_height - total_height) / 2.0
    cursor = lower
    for artist, height in zip(artists, heights):
        artist.set_position((0.075, (cursor + height / 2.0) / figure_height))
        cursor += height + gap


def single_panel(path: Path, csv_path: Path, mesh_label: str, color: str) -> None:
    data = load_csv(csv_path)
    x, healthy_curve, fault_curve, _, _ = local_window(data)
    fig, ax = plt.subplots(figsize=(14.0 / CM_PER_INCH, 8.0 / CM_PER_INCH), dpi=DPI)
    healthy_line, = ax.plot(x, healthy_curve / 1e8, color="#000000", linewidth=0.90, label="健康")
    fault_line, = ax.plot(x, fault_curve / 1e8, color=color, linewidth=1.45, label="故障")
    ax.axvline(CENTER_CYCLE, color="#8e8e8e", linestyle="--", linewidth=0.60, zorder=0)
    ax.set_xlim(0, 6)
    ax.set_xticks(range(7))
    ax.set_xlabel("啮合相位（周期）", fontproperties=SIMSUN, labelpad=5)
    # Curves are already plotted as y / 1e8.  State that scale only once in
    # the y-axis title and suppress Matplotlib's own offset/scientific text.
    add_mixed_ylabel(fig, mesh_label)
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax.yaxis.get_offset_text().set_visible(False)
    ax.tick_params(axis="both", labelsize=11, width=0.75, length=4.0, pad=3.0)
    for text in (*ax.get_xticklabels(), *ax.get_yticklabels()):
        text.set_fontproperties(TIMES)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    # The legend lives in dedicated top whitespace rather than over the data.
    fig.legend([healthy_line, fault_line], ["健康", "故障"], ncol=2, loc="upper center",
               bbox_to_anchor=(0.53, 0.99), prop=SIMSUN, frameon=False,
               handlelength=2.0, handletextpad=0.45, columnspacing=1.2)
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.19, top=0.84)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    configure_style()
    mpl.rcParams.update({
        "mathtext.fontset": "custom",
        "mathtext.rm": "Times New Roman",
        "mathtext.it": "Times New Roman:italic",
        "mathtext.bf": "Times New Roman:bold",
    })
    single_panel(run_dir / "Fig2_10a_planet_fault_SP_TVMS_local_14cmx8cm.png",
                 run_dir / "planet_fault_sp_p1_q050.csv", "SP", "#1f4e79")
    single_panel(run_dir / "Fig2_10b_planet_fault_PR_TVMS_local_14cmx8cm.png",
                 run_dir / "planet_fault_pr_p1_q050.csv", "PR", "#a64b00")
    print("Exported separate 14 cm x 8 cm thesis figures with 11-point text.")


if __name__ == "__main__":
    main()
