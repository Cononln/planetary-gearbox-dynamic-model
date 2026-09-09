"""Engineering schematic of the root-crack geometry used by the TVMS model.

The figure follows the potential-energy single-tooth coordinate convention.  It
is an original rendering of the analytical variables, not a copied source image.
"""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "results" / "tvms" / ".matplotlib_crack_geometry"))

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import FancyArrowPatch, Arc


OUT = ROOT / "results" / "tvms"
CM_PER_INCH = 2.54
TIMES = FontProperties(fname=r"C:\Windows\Fonts\times.ttf", size=11)
TIMES_SMALL = FontProperties(fname=r"C:\Windows\Fonts\times.ttf", size=9.5)
SIMSUN_SMALL = FontProperties(fname=r"C:\Windows\Fonts\simsun.ttc", size=9.5)
BLACK, GUIDE, RED = "#1e2329", "#626c78", "#c63232"


def arr(ax: plt.Axes, a: tuple[float, float], b: tuple[float, float], **kw: object) -> None:
    opts = dict(arrowstyle="->", mutation_scale=8, linewidth=0.70, color=BLACK)
    opts.update(kw)
    ax.add_patch(FancyArrowPatch(a, b, **opts))


def dim(ax: plt.Axes, a: tuple[float, float], b: tuple[float, float], **kw: object) -> None:
    opts = dict(arrowstyle="<->", mutation_scale=7, linewidth=0.65, color=BLACK)
    opts.update(kw)
    ax.add_patch(FancyArrowPatch(a, b, **opts))


def guide(ax: plt.Axes, a: tuple[float, float], b: tuple[float, float], **kw: object) -> None:
    opts = dict(color=GUIDE, linewidth=0.66, linestyle=(0, (3.5, 2.5)), zorder=0)
    opts.update(kw)
    ax.plot([a[0], b[0]], [a[1], b[1]], **opts)


def create() -> plt.Figure:
    mpl.rcParams.update({
        "font.family": "Times New Roman", "mathtext.fontset": "custom",
        "mathtext.rm": "Times New Roman", "mathtext.it": "Times New Roman:italic",
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })
    fig, ax = plt.subplots(figsize=(14 / CM_PER_INCH, 8 / CM_PER_INCH), dpi=600)
    fig.subplots_adjust(left=0.02, right=0.985, bottom=0.04, top=0.97)
    ax.set(xlim=(0, 14), ylim=(0, 8)); ax.axis("off")

    # Key analytical points: O is gear centre, P is the local root reference,
    # C is the crack origin and M is the contact point.
    O, P, C, M = (1.05, 3.13), (6.18, 3.13), (6.93, 5.22), (11.25, 5.10)
    tooth_left_base, tooth_left_tip = (6.20, 3.23), (6.62, 5.45)
    tooth_right_base, tooth_right_tip = (8.08, 3.20), (8.57, 5.42)

    # Centreline and geometric construction rays are kept light and do not
    # compete with the actual tooth, crack, and force vectors.
    guide(ax, O, (12.45, 3.13), linewidth=0.72)
    for end in (tooth_left_base, tooth_left_tip, tooth_right_base, tooth_right_tip, M):
        guide(ax, O, end, linewidth=0.62)
    ax.text(O[0] - 0.30, O[1] - 0.22, r"$O$", fontproperties=TIMES)
    ax.text(5.97, 2.73, r"$r_b$", fontproperties=TIMES_SMALL, color="#4e5965")
    ax.text(4.50, 3.35, r"$r_f$", fontproperties=TIMES_SMALL, color="#4e5965")
    ax.text(5.00, 3.65, r"$\alpha_1$", fontproperties=TIMES_SMALL)
    ax.text(5.52, 2.68, r"$\alpha_F$", fontproperties=TIMES_SMALL)
    ax.add_patch(Arc(O, 3.15, 3.15, theta1=0, theta2=23, linewidth=0.58, color=GUIDE))
    ax.add_patch(Arc(O, 3.95, 3.95, theta1=-16, theta2=0, linewidth=0.58, color=GUIDE))

    # Simplified involute flanks and top land; enough detail to support the
    # potential-energy coordinates without implying unmeasured fillet geometry.
    ax.plot([tooth_left_base[0], 6.42, tooth_left_tip[0]],
            [tooth_left_base[1], 4.47, tooth_left_tip[1]], color=BLACK, linewidth=1.05)
    ax.plot([tooth_right_base[0], 8.38, tooth_right_tip[0]],
            [tooth_right_base[1], 4.50, tooth_right_tip[1]], color=BLACK, linewidth=1.05)
    ax.plot([tooth_left_tip[0], tooth_right_tip[0]], [tooth_left_tip[1], tooth_right_tip[1]],
            color=BLACK, linewidth=1.05)
    ax.plot([tooth_left_base[0], tooth_right_base[0]], [tooth_left_base[1], tooth_right_base[1]],
            color=BLACK, linewidth=0.85)

    # Red straight root crack.  Its extension and inclination form the two
    # fault parameters q_c and gamma used by the implemented crack model.
    crack_end = (5.76, 7.08)
    ax.plot([C[0], crack_end[0]], [C[1], crack_end[1]], color=RED, linewidth=1.28,
            linestyle=(0, (5, 2.4)), zorder=4)
    ax.plot([C[0] - 0.12, crack_end[0] - 0.12], [C[1] + 0.04, crack_end[1] + 0.04],
            color=RED, linewidth=0.58, linestyle=(0, (2.2, 2.0)), zorder=3)
    dim(ax, (6.80, 5.40), (5.66, 7.16), color=RED, linewidth=0.72)
    ax.text(5.95, 6.41, r"$q_c$", fontproperties=TIMES, color=RED)
    ax.add_patch(Arc(C, 0.85, 0.85, theta1=78, theta2=123, linewidth=0.65, color=RED))
    ax.text(6.24, 5.91, r"$\gamma$", fontproperties=TIMES_SMALL, color=RED)
    ax.annotate("裂纹起点", xy=C, xytext=(3.52, 6.28), fontproperties=SIMSUN_SMALL, color="#4e5965",
                arrowprops=dict(arrowstyle="->", color=GUIDE, linewidth=0.58))

    # Contact force and its two components.
    arr(ax, M, (11.95, 5.82), linewidth=0.90)
    arr(ax, M, (12.02, 5.10), color=RED, linewidth=0.90)
    arr(ax, M, (11.25, 5.87), color=RED, linewidth=0.90)
    ax.text(12.03, 5.84, r"$F$", fontproperties=TIMES)
    ax.text(12.18, 5.05, r"$F_a$", fontproperties=TIMES_SMALL, color=RED)
    ax.text(11.35, 5.91, r"$F_b$", fontproperties=TIMES_SMALL, color=RED)

    # Horizontal and vertical dimensions defining the integration section.
    top_y, base_y = 6.42, 3.15
    guide(ax, (C[0], C[1]), (C[0], top_y), linewidth=0.50)
    guide(ax, (M[0], M[1]), (M[0], top_y), linewidth=0.50)
    dim(ax, (C[0], top_y), (M[0], top_y))
    ax.text(8.85, 6.60, r"$d$", fontproperties=TIMES_SMALL, ha="center")
    guide(ax, (7.56, 5.32), (7.56, top_y), linewidth=0.48)
    dim(ax, (C[0], 5.88), (7.56, 5.88))
    ax.text(7.28, 6.07, r"$x$", fontproperties=TIMES_SMALL)
    guide(ax, (M[0], base_y), M, linewidth=0.48)
    dim(ax, (M[0] + 0.43, base_y), (M[0] + 0.43, M[1]))
    ax.text(M[0] + 0.58, 4.10, r"$h$", fontproperties=TIMES_SMALL)

    # h_c, h_x and h_q are shown from the identical datum, avoiding ambiguous
    # arrows and matching the root-crack section definition in the model.
    for x, y, label in ((8.58, 4.13, r"$h_c$"), (7.84, 4.54, r"$h_x$"), (6.94, 3.82, r"$h_q$")):
        guide(ax, (x, base_y), (x, y), linewidth=0.46)
        dim(ax, (x, base_y), (x, y), linewidth=0.62)
        ax.text(x + 0.13, (base_y + y) / 2, label, fontproperties=TIMES_SMALL, va="center")
    ax.text(8.52, 2.57, "齿轮中心线", fontproperties=SIMSUN_SMALL, color="#4e5965")
    return fig


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig = create()
    stem = OUT / "Fig2_9_root_crack_geometry"
    fig.savefig(stem.with_suffix(".png"), dpi=600)
    fig.savefig(stem.with_suffix(".svg"))
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)
    print(f"Exported {stem}.png/.svg/.pdf")


if __name__ == "__main__":
    main()
