#!/usr/bin/env python3
"""中文证据图：齿圈相对齿序、故障增量幅值与模态峰值先后。

本脚本只读取已经保存的事件指标和四测点时域记录，不重算 18DOF、TVMS 或
ANSYS 模态路径。齿号是模型相对编号：约定仿真 t=0 时，P1 位于 SENSOR_0
径向线上，R1 开始进入 PR 啮合，齿圈顺时针编号为 R1 -> R2 -> ... -> R84。
ANSYS 资产没有逐齿实物编号，因此图中不会把模型编号冒充 CAD 齿号。
"""

from __future__ import annotations

import sys
from collections import OrderedDict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Circle

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))
from figure_qa import require_matplotlib_panel_alignment


RESULTS = ROOT / "results" / "fe_causal_q090_corrected_registration_20260928"
EVENT_CSV = RESULTS / "corrected_event_metrics.csv"
SIGNAL_CSV = RESULTS / "corrected_registration_time_signals.csv"
OUT = RESULTS / "mapping_evidence"
TOOTH_CSV = OUT / "model_relative_ring_tooth_sequence.csv"

N_RING = 84
SENSOR_ORDER = [0.0, 90.0, 120.0, 240.0]
SENSOR_COLORS = {
    0.0: "#1769aa",
    90.0: "#d95f02",
    120.0: "#2a9d55",
    240.0: "#7b4ab5",
}
STATE_ORDER = OrderedDict(
    [(1, ("E1", 177.2241209)), (2, ("E2", 267.2259520)),
     (3, ("E3", 357.2222899)), (4, ("E4", 87.2241209))]
)
STATE_COLORS = {1: "#4c78a8", 2: "#f58518", 3: "#54a24b", 4: "#e45756"}

FONT_CN_PATH = Path(r"C:/Windows/Fonts/simsun.ttc")
FONT_EN_PATH = Path(r"C:/Windows/Fonts/times.ttf")
FONT_CN = FontProperties(fname=str(FONT_CN_PATH)) if FONT_CN_PATH.exists() else FontProperties()
FONT_EN = FontProperties(fname=str(FONT_EN_PATH)) if FONT_EN_PATH.exists() else FontProperties(family="Times New Roman")

mpl.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "SimSun", "DejaVu Serif"],
        "font.sans-serif": ["SimSun"],
        "font.size": 8.0,
        "axes.titlesize": 9.0,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.2,
        "ytick.labelsize": 7.2,
        "legend.fontsize": 7.2,
        "axes.linewidth": 0.75,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.unicode_minus": False,
    }
)


def apply_cn_text(obj) -> None:
    """Use SimSun for Chinese text while retaining editable PDF/SVG text."""

    try:
        obj.set_fontproperties(FONT_CN)
    except AttributeError:
        pass


def save_figure(fig: mpl.figure.Figure, stem: Path, multi_panel: bool = False) -> None:
    fig.canvas.draw()
    if multi_panel:
        require_matplotlib_panel_alignment(
            fig,
            json_out=f"{stem}.alignment.json",
            overlay_svg=f"{stem}.alignment.svg",
            tolerance_pt=1.5,
            gutter_tolerance_pt=1.5,
            strict=True,
        )
    fig.savefig(f"{stem}.svg", bbox_inches="tight")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    events = pd.read_csv(EVENT_CSV)
    signals = pd.read_csv(SIGNAL_CSV)
    teeth = pd.read_csv(TOOTH_CSV)
    events["event_number"] = events["event_number"].astype(int)
    events["state"] = (events["event_number"] - 1) % 4 + 1
    numeric = [
        "active_planet", "contact_angle_global_deg", "sensor_angle_local_deg",
        "sensor_angle_global_deg", "contact_sensor_angular_distance_deg",
        "source_event_start_time_s", "source_event_peak_time_s",
        "source_event_end_time_s", "first_5pct_peak_crossing_lag_s",
        "five_pct_energy_lag_s", "peak_lag_s", "increment_peak_abs_ms2",
    ]
    for col in numeric:
        events[col] = events[col].astype(float)
    if len(teeth) != 48 or teeth["event_number"].tolist() != list(range(1, 49)):
        raise ValueError("The MAT-derived ring-tooth table must contain events 1..48")
    return events, signals, teeth


def nominal_fe_angle_for_ring_tooth(tooth: int) -> float:
    return float((90.0 - (tooth - 1) * 360.0 / N_RING) % 360.0)


def event_summary(events: pd.DataFrame, teeth: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for event_number, group in events.groupby("event_number", sort=True):
        group = group.sort_values("sensor_angle_local_deg")
        first = group.iloc[0]
        nearest = group.loc[group["contact_sensor_angular_distance_deg"].idxmin()]
        strongest = group.loc[group["increment_peak_abs_ms2"].idxmax()]
        s0 = group.loc[group["sensor_angle_local_deg"].eq(0.0)].iloc[0]
        s90 = group.loc[group["sensor_angle_local_deg"].eq(90.0)].iloc[0]
        tooth = teeth.loc[teeth["event_number"].eq(event_number)].iloc[0]
        if int(tooth["active_planet"]) != int(first["active_planet"]):
            raise ValueError(f"Event {event_number}: MAT tooth table and metric planet differ")
        if not np.isclose(float(tooth["source_event_peak_time_s"]),
                          float(first["source_event_peak_time_s"]), rtol=0, atol=1e-10):
            raise ValueError(f"Event {event_number}: MAT tooth table and metric time differ")
        previous_label, entry_label = str(tooth["peak_ring_tooth_pair"]).split("/")
        previous, entry = int(previous_label[1:]), int(entry_label[1:])
        # Positive means SENSOR_0 reaches its event-window modal peak earlier.
        s0_lead_ms = (float(s90["peak_lag_s"]) - float(s0["peak_lag_s"])) * 1000.0
        rows.append(
            {
                "event_number": int(event_number),
                "state": int(first["state"]),
                "active_planet": int(first["active_planet"]),
                "contact_angle_global_deg": float(first["contact_angle_global_deg"]),
                "source_event_start_time_s": float(first["source_event_start_time_s"]),
                "source_event_peak_time_s": float(first["source_event_peak_time_s"]),
                "source_event_end_time_s": float(first["source_event_end_time_s"]),
                "nearest_sensor_local_deg": float(nearest["sensor_angle_local_deg"]),
                "nearest_distance_deg": float(nearest["contact_sensor_angular_distance_deg"]),
                "strongest_sensor_local_deg": float(strongest["sensor_angle_local_deg"]),
                "strongest_peak_abs_ms2": float(strongest["increment_peak_abs_ms2"]),
                "nearest_equals_strongest": bool(
                    nearest["sensor_angle_local_deg"] == strongest["sensor_angle_local_deg"]
                ),
                "S0_peak_abs_ms2": float(s0["increment_peak_abs_ms2"]),
                "S90_peak_abs_ms2": float(s90["increment_peak_abs_ms2"]),
                "S0_peak_lag_s": float(s0["peak_lag_s"]),
                "S90_peak_lag_s": float(s90["peak_lag_s"]),
                "S0_peak_lead_ms": s0_lead_ms,
                "S0_crossing_lead_ms": (
                    float(s90["first_5pct_peak_crossing_lag_s"])
                    - float(s0["first_5pct_peak_crossing_lag_s"])
                ) * 1000.0,
                "peak_previous_ring_tooth": previous,
                "peak_new_ring_tooth": entry,
                "peak_ring_tooth_pair": str(tooth["peak_ring_tooth_pair"]),
                "whole_event_ring_tooth_set": str(tooth["whole_event_ring_tooth_set"]),
                "peak_new_tooth_nominal_fe_deg": nominal_fe_angle_for_ring_tooth(entry),
            }
        )
    return pd.DataFrame(rows)


def write_sequence_csv(summary: pd.DataFrame) -> Path:
    out = OUT / "model_relative_ring_tooth_sequence_cn.csv"
    table = summary.copy()
    table["齿序约定"] = "t=0：P1/SENSOR_0径向线，R1开始进入PR；顺时针 R1→…→R84"
    table = table.rename(
        columns={
            "event_number": "事件编号",
            "active_planet": "活跃行星",
            "source_event_start_time_s": "事件起点_保存段_s",
            "source_event_peak_time_s": "事件峰值_保存段_s",
            "source_event_end_time_s": "事件终点_保存段_s",
            "contact_angle_global_deg": "FE接触角_deg",
            "peak_previous_ring_tooth": "峰值时刻上一枚齿",
            "peak_new_ring_tooth": "峰值时刻新入啮齿",
            "peak_ring_tooth_pair": "峰值时刻齿圈齿对",
            "whole_event_ring_tooth_set": "事件全过程齿圈齿",
            "nearest_sensor_local_deg": "最近测点_local_deg",
            "strongest_sensor_local_deg": "最大峰值测点_local_deg",
            "nearest_distance_deg": "最近角距离_deg",
            "strongest_peak_abs_ms2": "四测点最大峰值_ms2",
            "S0_peak_abs_ms2": "S0峰值_ms2",
            "S90_peak_abs_ms2": "S90峰值_ms2",
            "S0_peak_lead_ms": "S0相对S90峰值领先_ms",
            "peak_new_tooth_nominal_fe_deg": "新入啮齿名义FE方位_deg",
        }
    )
    table.to_csv(out, index=False, encoding="utf-8-sig")
    return out


def plot_geometry(events: pd.DataFrame, summary: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.7))
    ax.set_aspect("equal")
    ax.set_xlim(-1.95, 1.95)
    ax.set_ylim(-1.62, 1.70)
    ax.axis("off")
    ax.add_patch(Circle((0, 0), 1.0, fill=False, lw=1.25, edgecolor="#2b2b2b"))
    # 84 个相对齿位：顺时针编号；只标出序列表中关键齿，避免密集文字覆盖。
    for tooth in range(1, N_RING + 1):
        theta = np.deg2rad(nominal_fe_angle_for_ring_tooth(tooth))
        r0, r1 = 1.005, 1.075 if tooth % 3 == 0 else 1.045
        ax.plot([r0 * np.cos(theta), r1 * np.cos(theta)],
                [r0 * np.sin(theta), r1 * np.sin(theta)], color="#8a8a8a", lw=0.30)

    sensor_global = {
        float(sensor): float(group.iloc[0]["sensor_angle_global_deg"])
        for sensor, group in events.groupby("sensor_angle_local_deg")
    }
    sensor_label_positions = {
        0.0: (0.0, 1.42, "center"), 90.0: (1.34, 0.05, "left"),
        120.0: (1.34, -0.78, "left"), 240.0: (-1.34, -0.78, "right"),
    }
    for sensor in SENSOR_ORDER:
        theta = np.deg2rad(sensor_global[sensor])
        x, y = 1.12 * np.cos(theta), 1.12 * np.sin(theta)
        ax.plot(x, y, marker="s", ms=6.0, color=SENSOR_COLORS[sensor], mec="white", mew=0.7)
        tx, ty, ha = sensor_label_positions[sensor]
        text = ax.text(tx, ty, f"测点 {int(sensor)}°\n(FE {sensor_global[sensor]:.1f}°)",
                       color=SENSOR_COLORS[sensor], fontsize=8.0, ha=ha, va="center")
        apply_cn_text(text)

    # The first four recurrent source events define the tooth sequence.
    for state, (label, _) in STATE_ORDER.items():
        row = summary.loc[summary["state"].eq(state)].iloc[0]
        contact = float(row["contact_angle_global_deg"])
        theta = np.deg2rad(contact)
        x, y = np.cos(theta), np.sin(theta)
        ax.plot([0.89 * x, 0.96 * x], [0.89 * y, 0.96 * y],
                color="#c62828", lw=0.7, ls="--", alpha=0.70)
        ax.plot(0.96 * x, 0.96 * y, marker="*", ms=11.0, color="#c62828", mec="white", mew=0.6)
        pair = str(row["peak_ring_tooth_pair"])
        ann = ax.text(0.70 * x, 0.70 * y, f"{label}\n{pair}\n{contact:.1f}°",
                      ha="center", va="center", fontsize=7.2, color="#a51e1e",
                      bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 1.2})
        apply_cn_text(ann)

    ax.annotate("顺时针齿序", xy=(-0.58, 0.82), xytext=(-0.98, 1.17),
                arrowprops={"arrowstyle": "->", "lw": 0.7, "color": "#333333"},
                fontsize=8.0, color="#333333")
    apply_cn_text(ax.texts[-1])
    title = ax.set_title("齿圈模型相对齿序与四种重复啮合位置", fontsize=10.2, pad=10)
    apply_cn_text(title)
    note1 = fig.text(0.5, 0.055,
                     "约定：仿真 t=0 时 P1 位于测点0°径向线上，R1 开始进入 PR 啮合；齿圈按顺时针 R1→R2→…→R84。",
                     ha="center", va="bottom", fontsize=7.2)
    note2 = fig.text(0.5, 0.028,
                     "红色星号为活跃行星的 PR 接触位置；标注为峰值时刻并联齿对。编号是模型相对编号，不是 ANSYS 实物刻号。",
                     ha="center", va="bottom", fontsize=7.2)
    apply_cn_text(note1); apply_cn_text(note2)
    fig.subplots_adjust(top=0.91, bottom=0.11)
    save_figure(fig, OUT / "Fig_FE_path_event_geometry_cn")


def plot_sequence_table(summary: pd.DataFrame) -> None:
    first = summary.head(12).copy()
    rows = []
    for _, row in first.iterrows():
        rows.append([
            f"{int(row['event_number'])}",
            f"P{int(row['active_planet'])}",
            f"{1000*row['source_event_peak_time_s']:.3f}",
            f"R{int(row['peak_previous_ring_tooth'])}",
            f"R{int(row['peak_new_ring_tooth'])}",
            str(row["whole_event_ring_tooth_set"]).replace(";", ","),
            f"{row['contact_angle_global_deg']:.2f}",
            f"S{int(row['nearest_sensor_local_deg'])}/S{int(row['strongest_sensor_local_deg'])}",
        ])
    fig, ax = plt.subplots(figsize=(9.4, 4.85))
    ax.axis("off")
    columns = ["事件", "行星", "峰值/ms", "上一并联齿", "新入啮齿",
               "全过程齿", "FE角/°", "近/强测点"]
    table = ax.table(cellText=rows, colLabels=columns, cellLoc="center", colLoc="center",
                     loc="upper center", bbox=[0.02, 0.17, 0.96, 0.72],
                     colWidths=[0.07, 0.08, 0.13, 0.105, 0.105, 0.24, 0.11, 0.16])
    table.auto_set_font_size(False)
    table.set_fontsize(7.6)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#a7a7a7")
        cell.set_linewidth(0.45)
        cell.get_text().set_fontproperties(FONT_CN)
        if r == 0:
            cell.set_facecolor("#e6eef5")
            cell.get_text().set_weight("bold")
        elif r % 4 == 0:
            cell.set_facecolor("#f6f6f6")
    title = ax.set_title("齿圈相对齿序表：保存段前 12 个太阳轮裂纹事件", fontsize=10.5, pad=13)
    apply_cn_text(title)
    note = fig.text(0.5, 0.085,
                    "保存段事件1：P2；峰值瞬间 R64 新入啮、R63 为上一并联齿，事件全过程还涉及 R65。",
                    ha="center", va="bottom", fontsize=7.7)
    note2 = fig.text(0.5, 0.045,
                     "齿号为模型相对编号；表中时间已舍弃前 0.2 s 过渡段。48 个事件的完整齿序见 CSV。",
                     ha="center", va="bottom", fontsize=7.7)
    apply_cn_text(note); apply_cn_text(note2)
    fig.subplots_adjust(top=0.88, bottom=0.14)
    save_figure(fig, OUT / "Fig_FE_path_event_sequence_cn")


def plot_amplitude_phase(summary: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.4, 6.0), sharex=True)
    x = summary["event_number"].to_numpy(int)
    a0 = 1000.0 * summary["S0_peak_abs_ms2"].to_numpy(float)
    a90 = 1000.0 * summary["S90_peak_abs_ms2"].to_numpy(float)
    ax = axes[0]
    ax.plot(x, a0, "o-", ms=2.6, lw=0.85, color=SENSOR_COLORS[0.0], label="测点0°")
    ax.plot(x, a90, "o-", ms=2.6, lw=0.85, color=SENSOR_COLORS[90.0], label="测点90°")
    ax.set_ylabel(r"故障增量峰值 |Δa| / 10$^{-3}$ m·s$^{-2}$", fontproperties=FONT_CN)
    ax.grid(alpha=0.18, lw=0.4)
    ax.set_title("(a) 0°/90°测点：逐事件故障增量峰值", fontproperties=FONT_CN, loc="left", pad=5)

    ax = axes[1]
    lead = summary["S0_peak_lead_ms"].to_numpy(float)
    colors = np.where(lead >= 0, SENSOR_COLORS[0.0], SENSOR_COLORS[90.0])
    ax.bar(x, lead, color=colors, width=0.72, alpha=0.82)
    ax.axhline(0, color="#333333", lw=0.65)
    ax.set_ylabel("S0相对S90峰值领先 / ms", fontproperties=FONT_CN)
    ax.set_xlabel("事件编号", fontproperties=FONT_CN)
    ax.set_title("(b) 峰值先后：正值表示测点0°先达峰，负值表示测点90°先达峰", fontproperties=FONT_CN, loc="left", pad=5)
    ax.grid(axis="y", alpha=0.18, lw=0.4)
    for ax in axes:
        ax.tick_params(direction="out", width=0.65, length=3)
        ax.set_xlim(0, 49)
    axes[-1].set_xticks(np.arange(0, 49, 4))
    fig.suptitle("齿圈传递路径的幅值与峰值先后证据", fontproperties=FONT_CN, fontsize=11, y=0.98)
    fig.legend(handles=axes[0].lines[:2], labels=["测点0°", "测点90°"],
               loc="upper center", bbox_to_anchor=(0.5, 0.935), ncol=2,
               frameon=False, prop=FONT_CN)
    note = fig.text(0.5, 0.025,
                    "48 个事件为同一确定性仿真的重复几何状态；峰值先后不是弹性波首达时间，幅值仅作相对比较。",
                    ha="center", va="bottom", fontsize=7.4)
    apply_cn_text(note)
    fig.subplots_adjust(left=0.13, right=0.985, bottom=0.13, top=0.84, hspace=0.29)
    save_figure(fig, OUT / "Fig_FE_path_amplitude_phase_cn", multi_panel=True)


def plot_statewise_amplitude_phase(summary: pd.DataFrame) -> None:
    """Statewise paired view that avoids implying a global 90-degree lead."""

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.65), constrained_layout=False)
    states = [1, 2, 3, 4]
    labels = []
    for state in states:
        first = summary.loc[summary["state"].eq(state)].iloc[0]
        # FE: sensor 0° is at global 90°; the user's angle grows clockwise.
        user_angle = (90.0 - float(first["contact_angle_global_deg"])) % 360.0
        labels.append(f"状态{state}\n≈{user_angle:.1f}°")

    # (a) Paired absolute peak amplitudes, retaining the event pairing.
    ax = axes[0]
    for i, state in enumerate(states):
        g = summary.loc[summary["state"].eq(state)].sort_values("event_number")
        jitter = np.linspace(-0.22, 0.22, len(g))
        a0 = 1000.0 * g["S0_peak_abs_ms2"].to_numpy(float)
        a90 = 1000.0 * g["S90_peak_abs_ms2"].to_numpy(float)
        for xj, y0, y90 in zip(jitter, a0, a90):
            ax.plot([i + xj - 0.07, i + xj + 0.07], [y0, y90],
                    color="#bdbdbd", lw=0.55, zorder=1)
        ax.scatter(i + jitter - 0.07, a0, s=18, color=SENSOR_COLORS[0.0], zorder=3,
                   label="测点0°" if i == 0 else None)
        ax.scatter(i + jitter + 0.07, a90, s=18, color=SENSOR_COLORS[90.0], zorder=3,
                   label="测点90°" if i == 0 else None)
    ax.set_xticks(range(4), labels, fontproperties=FONT_CN)
    ax.set_ylabel(r"故障增量峰值 |Δa| / 10$^{-3}$ m·s$^{-2}$", fontproperties=FONT_CN)
    ax.set_title("(a) 按啮合状态配对比较幅值", fontproperties=FONT_CN, loc="left", pad=5)
    ax.grid(axis="y", alpha=0.18, lw=0.4)

    # (b) Signed peak-time difference. Negative means 90° peak is earlier.
    ax = axes[1]
    for i, state in enumerate(states):
        g = summary.loc[summary["state"].eq(state)].sort_values("event_number")
        jitter = np.linspace(-0.22, 0.22, len(g))
        delta = 1000.0 * (g["S90_peak_lag_s"] - g["S0_peak_lag_s"]).to_numpy(float)
        ax.scatter(i + jitter, delta, s=18, color=STATE_COLORS[state], alpha=0.82,
                   edgecolor="white", linewidth=0.25, zorder=3)
        ax.plot([i - 0.25, i + 0.25], [np.median(delta), np.median(delta)],
                color="#222222", lw=1.1, zorder=4)
    ax.axhline(0, color="#333333", lw=0.65)
    ax.set_xticks(range(4), labels, fontproperties=FONT_CN)
    ax.set_ylabel(r"Δt$_{peak}$ = t$_{90}$ - t$_{0}$ / ms", fontproperties=FONT_CN)
    ax.set_title("(b) 峰值时间差：负值表示90°先出现峰值", fontproperties=FONT_CN, loc="left", pad=5)
    ax.grid(axis="y", alpha=0.18, lw=0.4)

    for ax in axes:
        ax.tick_params(direction="out", width=0.65, length=3)
    fig.suptitle("按四种实际重复接触相位展示 0°/90°测点响应", fontproperties=FONT_CN, fontsize=11, y=0.98)
    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, frameon=False, prop=FONT_CN, ncol=2,
               loc="upper left", bbox_to_anchor=(0.095, 0.91))
    note = fig.text(0.5, 0.018,
                    "状态每4个事件重复一次，并非人为等宽分箱；角度按用户坐标（0°向上、顺时针为正）。",
                    ha="center", va="bottom", fontsize=7.3)
    apply_cn_text(note)
    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.22, top=0.81,
                        wspace=0.27)
    save_figure(fig, OUT / "Fig_FE_path_statewise_amplitude_phase_cn", multi_panel=True)


def plot_phase_lead_cn(events: pd.DataFrame, signals: pd.DataFrame, summary: pd.DataFrame) -> None:
    time = signals["time_s"].to_numpy(float)
    y0 = signals["q0900um_minus_healthy_sensor_000deg_acceleration_ms2"].to_numpy(float) * 1000.0
    y90 = signals["q0900um_minus_healthy_sensor_090deg_acceleration_ms2"].to_numpy(float) * 1000.0
    selected = [4, 8]
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.0), sharex=True, sharey=True)
    ymax = 0.0
    records = []
    for ev in selected:
        rows = events.loc[events["event_number"].eq(ev)]
        start = float(rows.iloc[0]["source_event_start_time_s"])
        keep = (time >= start) & (time <= start + 0.008)
        t_ms = (time[keep] - start) * 1000.0
        yy0, yy90 = y0[keep], y90[keep]
        ymax = max(ymax, float(np.max(np.abs(np.r_[yy0, yy90]))))
        records.append((t_ms, yy0, yy90, rows))
    for ax, ev, (t_ms, yy0, yy90, rows) in zip(axes, selected, records):
        r0 = rows.loc[rows["sensor_angle_local_deg"].eq(0.0)].iloc[0]
        r90 = rows.loc[rows["sensor_angle_local_deg"].eq(90.0)].iloc[0]
        ax.plot(t_ms, yy0, color=SENSOR_COLORS[0.0], lw=0.9, label="测点0°")
        ax.plot(t_ms, yy90, color=SENSOR_COLORS[90.0], lw=0.9, label="测点90°")
        p0, p90 = float(r0["peak_lag_s"]) * 1000.0, float(r90["peak_lag_s"]) * 1000.0
        c0 = float(r0["first_5pct_peak_crossing_lag_s"]) * 1000.0
        c90 = float(r90["first_5pct_peak_crossing_lag_s"]) * 1000.0
        # Dashed lines mark first threshold crossing; dotted lines mark the
        # event-window maximum absolute modal-ringdown peak.
        ax.axvline(c0, color=SENSOR_COLORS[0.0], ls="--", lw=0.65, alpha=0.55)
        ax.axvline(c90, color=SENSOR_COLORS[90.0], ls="--", lw=0.65, alpha=0.55)
        ax.axvline(p0, color=SENSOR_COLORS[0.0], ls=":", lw=0.9)
        ax.axvline(p90, color=SENSOR_COLORS[90.0], ls=":", lw=0.9)
        i0, i90 = int(np.argmin(np.abs(t_ms - p0))), int(np.argmin(np.abs(t_ms - p90)))
        ax.scatter([t_ms[i0], t_ms[i90]], [yy0[i0], yy90[i90]],
                   c=[SENSOR_COLORS[0.0], SENSOR_COLORS[90.0]], s=18, zorder=6)
        ax.axhline(0.0, color="#777777", lw=0.45)
        row = summary.loc[summary["event_number"].eq(ev)].iloc[0]
        lead = float(row["S0_peak_lead_ms"])
        order = "S0先达峰" if lead > 0 else "S90先达峰"
        title = ax.set_title(
            f"事件{ev}（{row['peak_ring_tooth_pair']}）：{order} {abs(lead):.3f} ms；"
            f"最大峰 S0/S90 = {1000*row['S0_peak_abs_ms2']:.3f}/{1000*row['S90_peak_abs_ms2']:.3f}\n"
            f"首次5%阈值 t0/t90 = {c0:.3f}/{c90:.3f} ms",
            fontproperties=FONT_CN, fontsize=8.4, loc="left", pad=6)
        ax.set_ylabel(r"故障增量 / 10$^{-3}$ m·s$^{-2}$", fontproperties=FONT_CN)
        ax.set_xlim(0.0, 8.0)
        ax.set_ylim(-1.25 * ymax, 1.25 * ymax)
        ax.grid(alpha=0.17, lw=0.4)
    axes[-1].set_xlabel("距源事件起点 / ms", fontproperties=FONT_CN)
    fig.suptitle("0°/90°测点局部故障波形与模态峰值先后", fontproperties=FONT_CN, fontsize=10.6, y=0.98)
    fig.legend(handles=axes[0].lines[:2], labels=["测点0°", "测点90°"], loc="upper center",
               bbox_to_anchor=(0.5, 0.935), ncol=2, frameon=False, prop=FONT_CN)
    note = fig.text(0.5, 0.018,
                    "浅虚线=首次超过5%阈值；点虚线=事件窗最大绝对峰。最大峰先后不等于有限波速首达时间。",
                    ha="center", va="bottom", fontsize=7.0)
    apply_cn_text(note)
    fig.subplots_adjust(top=0.81, bottom=0.11, left=0.12, right=0.985, hspace=0.32)
    save_figure(fig, OUT / "Fig_FE_path_phase_lead_S0_S90_cn", multi_panel=True)


def write_notes(summary: pd.DataFrame) -> None:
    path = OUT / "中文说明_齿序_幅值_相位.txt"
    nearest_count = int(summary["nearest_equals_strongest"].sum())
    with path.open("w", encoding="utf-8") as f:
        f.write("齿圈传递路径中文证据图说明\n")
        f.write("================================\n")
        f.write("1. 齿号约定\n")
        f.write("仿真 t=0 时，P1 位于 SENSOR_0 径向线上，R1 开始进入 PR 啮合；齿圈按顺时针 R1→R2→…→R84。\n")
        f.write("ANSYS 模态资产没有逐齿实物编号，因此 R1~R84 是模型相对齿号，不是 CAD 刻号。\n")
        f.write("保存段已舍弃原仿真的 0.2 s；表中事件1是保存段的第一个完整太阳轮裂纹事件。\n")
        f.write("事件1峰值齿对为 R63/R64，事件全过程涉及 R63/R64/R65；前四事件峰值齿对依次为 R63/R64、R42/R43、R21/R22、R84/R1。\n")
        f.write("状态1~4不是将圆周人为划为四个90°扇区，而是将实际故障事件每4次重复出现的接触相位归组；四个相位依次相隔约90°。\n")
        f.write("将FE角度换算到用户约定（测点0°径向向上，顺时针为正）后，状态1~4接触位置约为272.8°、182.8°、92.8°、2.8°。\n")
        f.write("2. 幅值与峰值先后\n")
        f.write(f"最近测点同时为四测点最大峰值测点：{nearest_count}/{len(summary)} 个事件。\n")
        f.write("正的 S0_peak_lead_ms 表示 SENSOR_0 的事件窗峰值比 SENSOR_90 更早；该指标是模态振铃峰差，不是波前首达时间。\n")
        f.write("当前 FE 路径由模态状态因果积分构成，没有显式有限波速、距离延迟或经验距离衰减项。幅值和峰值先后还受振型节点/反节点、共振/反共振、多模态相消/相长、边界反射及三路 PR 力共同激励影响。\n")
        f.write("因此图中的结果支持‘角距离有统计趋势，但不能逐事件用距离唯一决定幅值和先后’，不支持简单的近点必强、近点必早。\n")
        f.write("当前 ANSYS 模态转 SI 后的广义质量归一化尚未核验，所以 m/s² 的绝对幅值不作为已标定定量结果，图中只比较同一计算设置下的相对幅值与峰值先后。\n")
        f.write("当前模型齿数 21/31/84 尚有标准同模数几何一致性问题，齿序图只能视为这套模型的相位记账，不能作为实物齿圈逐齿编号证据。\n")
        f.write("四状态分组图：同一事件的0°/90°故障增量绝对峰值配对；右图画 t90-t0 的逐事件数值和各组中位数，负值表示90°通道的最大模态峰先出现。\n")
        for state in (1, 2, 3, 4):
            g = summary.loc[summary["state"].eq(state)]
            n90_early = int((g["S90_peak_lag_s"] < g["S0_peak_lag_s"]).sum())
            n90_stronger = int((g["S90_peak_abs_ms2"] > g["S0_peak_abs_ms2"]).sum())
            f.write(f"状态{state}：90°峰值较早 {n90_early}/{len(g)}；90°幅值较大 {n90_stronger}/{len(g)}。\n")
        f.write("四组各12个事件来自同一确定性仿真反复经过的几何状态，不是12次独立试验。\n")
        f.write("3. 数据完整性\n")
        f.write("所有曲线直接读取保存的健康/故障差值和逐事件指标；为在图中使用10^-3刻度，仅做×1000单位换算。未平移、未平滑、未人为添加冲击/谱线。\n")
        f.write("完整48事件进入逐事件幅值和峰值先后图；齿序表展示前12事件以便阅读，但完整48行已另存CSV；局部时域图取事件4和8，分别展示接触点接近0°时最大峰由90°先出现与由0°先出现的两种情况，并叠加首次5%阈值。\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    events, signals, teeth = load_data()
    summary = event_summary(events, teeth)
    summary.to_csv(OUT / "event_geometry_phase_summary_cn.csv", index=False, encoding="utf-8-sig")
    write_sequence_csv(summary)
    plot_geometry(events, summary)
    plot_sequence_table(summary)
    plot_amplitude_phase(summary)
    plot_statewise_amplitude_phase(summary)
    plot_phase_lead_cn(events, signals, summary)
    write_notes(summary)
    print(f"中文证据图已写入: {OUT}")
    print(f"事件数={len(summary)}; 最近测点=最大峰值: {int(summary['nearest_equals_strongest'].sum())}/{len(summary)}")
    print("前四事件齿圈峰值齿对: " + " -> ".join(summary.head(4)["peak_ring_tooth_pair"].tolist()))


if __name__ == "__main__":
    main()
