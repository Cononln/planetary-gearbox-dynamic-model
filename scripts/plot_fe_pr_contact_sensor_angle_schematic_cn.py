#!/usr/bin/env python3
"""Plot the registered PR-contact and sensor azimuths as a schematic.

The four event positions come from the saved 0.90 mm sun-crack path diagnostic.
They denote distinct times; the diagram is not a simultaneous contact map or
an as-built dimensional drawing. No dynamics or FE data are recomputed here.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Circle, Wedge

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))
from figure_qa import require_matplotlib_panel_alignment


SOURCE = ROOT / "results" / "fe_causal_q090_corrected_registration_20260928"
EVENT_CSV = SOURCE / "threshold50_20260929" / "event_50pct_time_amplitude_0_90.csv"
SENSOR_CSV = SOURCE / "response_summary.csv"
OUTPUT = SOURCE / "angle_mapping_schematic_20260929"
STEM = OUTPUT / "Fig_FE_PR_contact_sensor_angle_schematic_cn"

CN_FILE = Path(r"C:/Windows/Fonts/simsun.ttc")
CN = FontProperties(fname=str(CN_FILE)) if CN_FILE.exists() else FontProperties(family="SimSun")
EN = FontProperties(family="Times New Roman")
mpl.rcParams.update({
    # Preferred Latin and Chinese faces; Arial is only a portability fallback.
    "font.family": ["Times New Roman", "SimSun", "Arial"],
    "font.size": 9.2,
    "axes.linewidth": 0.7,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

SENSOR_BLUE = "#276a9c"
CONTACT_RED = "#c94b43"
INK = "#28323b"
MUTED = "#6b7680"


def xy(angle_deg: float, radius: float) -> tuple[float, float]:
    """User convention: 0° points up, angle increases clockwise."""
    rad = np.deg2rad(angle_deg)
    return radius * float(np.sin(rad)), radius * float(np.cos(rad))


def source_angles() -> tuple[dict[int, float], dict[int, float]]:
    grouped: dict[int, list[float]] = {i: [] for i in range(1, 5)}
    with EVENT_CSV.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            grouped[int(row["state"])].append(float(row["contact_angle_user_deg"]))
    if [len(grouped[i]) for i in range(1, 5)] != [12] * 4:
        raise ValueError("Expected 12 saved events at each of four phases.")
    contacts = {i: float(np.median(grouped[i])) for i in range(1, 5)}

    sensors: dict[int, float] = {}
    with SENSOR_CSV.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["registration"] == "corrected" and row["case"] == "healthy":
                label = int(round(float(row["sensor_local_deg"])))
                sensors[label] = float(row["sensor_local_deg"])
    if set(sensors) != {0, 90, 120, 240}:
        raise ValueError("Saved FE asset must contain sensor labels 0°, 90°, 120°, 240°.")
    if not (abs(contacts[4]) < 5 and abs(contacts[3] - 90) < 5):
        raise ValueError("The saved corrected-angle contact registration changed.")
    return contacts, sensors


def make_figure(contacts: dict[int, float], sensors: dict[int, float]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.20, 5.15), facecolor="white")
    ax = fig.add_axes([0.035, 0.12, 0.565, 0.77])
    ax.set_aspect("equal")
    ax.set_xlim(-2.40, 2.40)
    ax.set_ylim(-2.35, 2.35)
    ax.axis("off")

    # Thin annulus: PR contact is on the inner face; sensor markers are outside.
    ax.add_patch(Wedge((0, 0), 1.63, 0, 360, width=0.20,
                       facecolor="#edf0f2", edgecolor=INK, linewidth=0.85))
    ax.add_patch(Circle((0, 0), 0.047, color=INK, zorder=3))
    ax.text(0, -0.23, "轮系中心", ha="center", va="center",
            fontproperties=CN, fontsize=8.5, color=MUTED)

    # The 0° radial line demonstrates aligned azimuth, not shared coordinates.
    x_start, y_start = xy(0, 1.53)
    x_end, y_end = xy(0, 1.88)
    ax.plot([x_start, x_end], [y_start, y_end], color=SENSOR_BLUE, lw=0.75,
            ls=(0, (3, 3)), alpha=0.65, zorder=1)
    ax.text(-1.73, 1.48, "齿圈内缘\nPR啮合位置", ha="right", va="center",
            fontproperties=CN, fontsize=8.3, color=MUTED)

    # Separate event markers intentionally superpose four DIFFERENT instants.
    for state, angle in contacts.items():
        x, y = xy(angle, 1.43)
        ax.add_patch(Circle((x, y), 0.075, facecolor=CONTACT_RED,
                            edgecolor="white", linewidth=0.8, zorder=5))
        tx, ty = xy(angle, 1.20)
        ax.text(tx, ty, str(state), color=CONTACT_RED, ha="center", va="center",
                fontproperties=EN, fontsize=8.8, weight="bold", zorder=6)

    for label, angle in sensors.items():
        x, y = xy(angle, 1.90)
        ax.scatter([x], [y], marker="s", s=72, c=SENSOR_BLUE,
                   edgecolors="white", linewidths=0.75, zorder=7)
    # Direct labels sit outside the sensors and do not overlap the ring.
    ax.text(0.0, 2.16, "S0  测点0°", fontproperties=CN, fontsize=9.2,
            color=SENSOR_BLUE, ha="center", va="bottom")
    ax.text(2.13, 0.00, "S90\n测点90°", fontproperties=CN, fontsize=8.5,
            color=SENSOR_BLUE, ha="left", va="center")
    ax.text(1.50, -1.34, "S120\n测点120°", fontproperties=CN, fontsize=8.3,
            color=SENSOR_BLUE, ha="left", va="top")
    ax.text(-1.50, -1.34, "S240\n测点240°", fontproperties=CN, fontsize=8.3,
            color=SENSOR_BLUE, ha="right", va="top")

    # Curved arrow conveys only the positive coordinate direction.
    ax.annotate("", xy=xy(53, 2.15), xytext=xy(25, 2.15),
                arrowprops={"arrowstyle": "-|>", "color": MUTED, "lw": 0.9,
                            "connectionstyle": "arc3,rad=-0.12"})

    fig.text(0.625, 0.855, "固定传感器", fontproperties=CN,
             fontsize=10.0, color=SENSOR_BLUE)
    fig.text(0.625, 0.815, "蓝色方块：0°、90°、120°、240°", fontproperties=CN,
             fontsize=8.5, color=INK)
    fig.text(0.625, 0.775, "顶部为0°，顺时针为正", fontproperties=CN,
             fontsize=8.5, color=MUTED)
    fig.text(0.625, 0.742, "不同时刻的PR接触位置", fontproperties=CN,
             fontsize=10.0, color=CONTACT_RED)
    for state, ypos in zip(range(1, 5), (0.694, 0.647, 0.600, 0.553)):
        fig.text(0.625, ypos, f"{state}  ≈ {contacts[state]:.1f}°",
                 fontproperties=CN, fontsize=9.0, color=INK)
    fig.text(0.625, 0.455, "关键区别", fontproperties=CN,
             fontsize=10.0, color=INK)
    fig.text(0.625, 0.405, "接触点4与测点S0近乎同方位，", fontproperties=CN,
             fontsize=8.6, color=INK)
    fig.text(0.625, 0.365, "但不在同一个物理位置。", fontproperties=CN,
             fontsize=8.6, color=INK)
    fig.text(0.625, 0.295, "180°与270°是接触方位；", fontproperties=CN,
             fontsize=8.6, color=INK)
    fig.text(0.625, 0.255, "当前没有同名传感器。", fontproperties=CN,
             fontsize=8.6, color=INK)

    fig.suptitle("齿圈PR接触方位与固定传感器位置示意",
                 fontproperties=CN, fontsize=11.3, y=0.965)
    fig.text(0.5, 0.055,
             "红点代表四次不同时间的故障事件，并非同时接触；方位按当前模型注册约定，几何尺寸未按比例绘制。",
             ha="center", va="center", fontproperties=CN, fontsize=7.8, color=MUTED)
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig, json_out=f"{STEM}.alignment.json", overlay_svg=f"{STEM}.alignment.svg",
        tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True,
    )
    fig.savefig(f"{STEM}.pdf", bbox_inches="tight")
    fig.savefig(f"{STEM}.svg", bbox_inches="tight")
    fig.savefig(f"{STEM}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{STEM}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    event_angles, sensor_angles = source_angles()
    make_figure(event_angles, sensor_angles)
    print(STEM)
