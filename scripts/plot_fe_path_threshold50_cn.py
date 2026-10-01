#!/usr/bin/env python3
"""Audit 0°/90° fault-increment half-peak times from saved FE path responses.

This is a post-processing diagnostic only.  It does not rerun the 18-DOF source
or FE modal path, and does not overwrite the historical 5% event metrics.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))
from figure_qa import require_matplotlib_panel_alignment


SOURCE = ROOT / "results" / "fe_causal_q090_corrected_registration_20260928"
EVENTS = SOURCE / "corrected_event_metrics.csv"
SIGNALS = SOURCE / "corrected_registration_time_signals.csv"
OUTPUT = SOURCE / "threshold50_20260929"
FRACTION = 0.50
PAD_S = 0.004
SENSOR_COLOR = {0: "#1769aa", 90: "#d95f02"}
STATE_COLOR = {1: "#4c78a8", 2: "#f58518", 3: "#54a24b", 4: "#e45756"}

CN_FONT_FILE = Path(r"C:/Windows/Fonts/simsun.ttc")
CN_FONT = FontProperties(fname=str(CN_FONT_FILE)) if CN_FONT_FILE.exists() else FontProperties()
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "SimSun", "DejaVu Serif"],
    "font.size": 8.0,
    "axes.titlesize": 9.0,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.2,
    "ytick.labelsize": 7.2,
    "axes.linewidth": 0.75,
    "mathtext.fontset": "stix",
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.unicode_minus": False,
})


def load_inputs() -> tuple[pd.DataFrame, np.ndarray, dict[int, np.ndarray]]:
    events = pd.read_csv(EVENTS)
    if events.shape[0] != 48 * 4:
        raise ValueError("Expected 48 events recorded at four sensors.")
    events = events.loc[events["sensor_angle_local_deg"].isin([0.0, 90.0])].copy()
    if events.groupby("event_number").size().ne(2).any():
        raise ValueError("Every event must have one 0° and one 90° row.")
    cols = [
        "time_s",
        "q0900um_minus_healthy_sensor_000deg_acceleration_ms2",
        "q0900um_minus_healthy_sensor_090deg_acceleration_ms2",
    ]
    signals = pd.read_csv(SIGNALS, usecols=cols)
    time = signals["time_s"].to_numpy(dtype=float)
    if not np.all(np.diff(time) > 0):
        raise ValueError("Signal time must increase strictly.")
    values = {
        0: signals[cols[1]].to_numpy(dtype=float),
        90: signals[cols[2]].to_numpy(dtype=float),
    }
    return events, time, values


def compute_half_peak_metrics(
    events: pd.DataFrame, time: np.ndarray, values: dict[int, np.ndarray]
) -> pd.DataFrame:
    rows: list[dict[str, float | int]] = []
    for number, pair in events.groupby("event_number", sort=True):
        first = pair.iloc[0]
        start = float(first["source_event_start_time_s"])
        end = float(first["source_event_end_time_s"])
        if not np.allclose(pair["source_event_start_time_s"], start, atol=1e-12):
            raise ValueError(f"Source start differs between channels for event {number}.")
        left = int(np.searchsorted(time, max(time[0], start - PAD_S), side="left"))
        right = int(np.searchsorted(time, min(time[-1], end + PAD_S), side="right"))
        if right <= left:
            raise ValueError(f"Empty event window for event {number}.")
        row: dict[str, float | int] = {
            "event_number": int(number),
            "state": (int(number) - 1) % 4 + 1,
            "contact_angle_user_deg": (90.0 - float(first["contact_angle_global_deg"])) % 360.0,
            "source_event_start_s": start,
            "source_event_end_s": end,
        }
        for angle in (0, 90):
            x = values[angle][left:right]
            t = time[left:right]
            peak = float(np.max(np.abs(x)))
            if peak <= 0.0:
                raise ValueError(f"Zero fault-increment peak for event {number}, sensor {angle}°.")
            stored = float(pair.loc[pair["sensor_angle_local_deg"].eq(angle),
                                    "increment_peak_abs_ms2"].iloc[0])
            if not np.isclose(peak, stored, rtol=2e-8, atol=1e-12):
                raise ValueError(f"Recomputed peak differs from source CSV at event {number}, sensor {angle}°.")
            crossing = np.flatnonzero((t >= start) & (np.abs(x) >= FRACTION * peak))
            if crossing.size == 0:
                raise ValueError(f"No half-peak crossing for event {number}, sensor {angle}°.")
            t_cross = float(t[int(crossing[0])])
            row[f"peak_abs_sensor_{angle}_ms2"] = peak
            row[f"threshold50_sensor_{angle}_ms2"] = FRACTION * peak
            row[f"first50_sensor_{angle}_time_s"] = t_cross
            row[f"first50_sensor_{angle}_lag_s"] = t_cross - start
        row["first50_t90_minus_t0_ms"] = 1000.0 * (
            float(row["first50_sensor_90_time_s"]) - float(row["first50_sensor_0_time_s"])
        )
        row["earlier_half_peak_sensor_deg"] = 0 if row["first50_t90_minus_t0_ms"] > 0 else 90
        row["larger_peak_sensor_deg"] = 0 if row["peak_abs_sensor_0_ms2"] > row["peak_abs_sensor_90_ms2"] else 90
        rows.append(row)
    return pd.DataFrame(rows)


def save_figure(fig: mpl.figure.Figure, stem: Path) -> None:
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig,
        json_out=f"{stem}.alignment.json",
        overlay_svg=f"{stem}.alignment.svg",
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        strict=True,
    )
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.svg", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def state_figure(metrics: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.45), constrained_layout=False)
    labels = [
        f"状态{i}\n≈{float(metrics.loc[metrics['state'].eq(i), 'contact_angle_user_deg'].median()):.1f}°"
        for i in range(1, 5)
    ]

    # Panel a is the timing result: negative values mean 90° reaches its own
    # half-peak earlier.  The common source event time cancels in the difference.
    ax = axes[0]
    for i in range(1, 5):
        group = metrics.loc[metrics["state"].eq(i)].sort_values("event_number")
        x = (i - 1) + np.linspace(-0.20, 0.20, len(group))
        y = group["first50_t90_minus_t0_ms"].to_numpy(float)
        ax.scatter(x, y, color=STATE_COLOR[i], s=19, alpha=0.84,
                   edgecolor="white", linewidth=0.25, zorder=3)
        ax.plot([i - 1.23, i - 0.77], [np.median(y)] * 2,
                color="#222222", lw=1.1, zorder=4)
    ax.axhline(0.0, color="#333333", lw=0.7)
    ax.set_xticks(range(4), labels, fontproperties=CN_FONT)
    ax.set_ylabel(r"50%峰值时间差：$t_{50,90}-t_{50,0}$ / ms", fontproperties=CN_FONT)
    ax.set_title("(a) 同一事件中谁先达到自身半峰值", fontproperties=CN_FONT, loc="left", pad=5)
    ax.grid(axis="y", alpha=0.18, lw=0.4)

    # Panel b preserves raw (non-normalized) fault-increment magnitude.
    ax = axes[1]
    for i in range(1, 5):
        group = metrics.loc[metrics["state"].eq(i)].sort_values("event_number")
        x = (i - 1) + np.linspace(-0.20, 0.20, len(group))
        y0 = 1000.0 * group["peak_abs_sensor_0_ms2"].to_numpy(float)
        y90 = 1000.0 * group["peak_abs_sensor_90_ms2"].to_numpy(float)
        for xx, a0, a90 in zip(x, y0, y90):
            ax.plot([xx - 0.055, xx + 0.055], [a0, a90], color="#bdbdbd", lw=0.55, zorder=1)
        ax.scatter(x - 0.055, y0, s=18, color=SENSOR_COLOR[0], zorder=3,
                   label="测点0°" if i == 1 else None)
        ax.scatter(x + 0.055, y90, s=18, color=SENSOR_COLOR[90], zorder=3,
                   label="测点90°" if i == 1 else None)
    ax.set_xticks(range(4), labels, fontproperties=CN_FONT)
    ax.set_ylabel(r"故障增量绝对峰值 / $(10^{-3}\,\mathrm{m\,s^{-2}})$", fontproperties=CN_FONT)
    ax.set_title("(b) 同一事件的故障增量幅值", fontproperties=CN_FONT, loc="left", pad=5)
    ax.grid(axis="y", alpha=0.18, lw=0.4)

    for ax in axes:
        ax.tick_params(direction="out", width=0.65, length=3)
    handles, legend_labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, legend_labels, frameon=False, prop=CN_FONT, ncol=2,
               loc="upper center", bbox_to_anchor=(0.5, 0.895))
    fig.suptitle("0°/90°故障响应：自身50%峰值时间与实际幅值", fontproperties=CN_FONT,
                 fontsize=10.7, y=0.98)
    fig.text(0.5, 0.018,
             "50%门槛由每个测点各自的故障增量峰值定义；负时间差表示90°先达半峰值，不能解释为波前首达。",
             ha="center", va="bottom", fontproperties=CN_FONT, fontsize=7.2)
    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.205, top=0.75, wspace=0.27)
    save_figure(fig, OUTPUT / "Fig_FE_path_statewise_50pct_amplitude_cn")


def local_figure(metrics: pd.DataFrame, time: np.ndarray, values: dict[int, np.ndarray]) -> None:
    selected = (4, 8)
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.1), sharex=True, sharey=True)
    bound = 0.0
    snippets = []
    for number in selected:
        row = metrics.loc[metrics["event_number"].eq(number)].iloc[0]
        start = float(row["source_event_start_s"])
        take = (time >= start) & (time <= start + 0.008)
        t_ms = 1000.0 * (time[take] - start)
        y0 = 1000.0 * values[0][take]
        y90 = 1000.0 * values[90][take]
        bound = max(bound, float(np.max(np.abs(np.r_[y0, y90]))))
        snippets.append((row, t_ms, y0, y90))
    for ax, number, (row, t_ms, y0, y90) in zip(axes, selected, snippets):
        ax.plot(t_ms, y0, color=SENSOR_COLOR[0], lw=0.85, label="测点0°")
        ax.plot(t_ms, y90, color=SENSOR_COLOR[90], lw=0.85, label="测点90°")
        for angle in (0, 90):
            cross_ms = 1000.0 * float(row[f"first50_sensor_{angle}_lag_s"])
            ax.axvline(cross_ms, color=SENSOR_COLOR[angle], lw=0.85, ls="--")
        ax.axhline(0.0, color="#777777", lw=0.45)
        dt_ms = float(row["first50_t90_minus_t0_ms"])
        earlier = "0°先达到" if dt_ms > 0 else "90°先达到"
        ax.set_title(
            f"事件{number}：{earlier}自身50%峰值，时间差 {abs(dt_ms):.3f} ms；"
            f"绝对峰值0°/90° = {1000*float(row['peak_abs_sensor_0_ms2']):.3f}/"
            f"{1000*float(row['peak_abs_sensor_90_ms2']):.3f}",
            fontproperties=CN_FONT, fontsize=8.2, loc="left", pad=6)
        ax.set_ylabel(r"故障增量 / $(10^{-3}\,\mathrm{m\,s^{-2}})$", fontproperties=CN_FONT)
        ax.set_xlim(0.0, 8.0)
        ax.set_ylim(-1.22 * bound, 1.22 * bound)
        ax.grid(alpha=0.16, lw=0.4)
    axes[-1].set_xlabel("距源事件起点 / ms", fontproperties=CN_FONT)
    fig.suptitle("接触位置靠近0°时的两次故障增量波形", fontproperties=CN_FONT,
                 fontsize=10.6, y=0.98)
    fig.legend(handles=axes[0].lines[:2], labels=["测点0°", "测点90°"],
               loc="upper center", bbox_to_anchor=(0.5, 0.935), ncol=2,
               frameon=False, prop=CN_FONT)
    fig.text(0.5, 0.017,
             "竖虚线=各测点首次达到自身事件窗最大绝对故障增量的50%；两图使用同一幅值尺度。",
             ha="center", va="bottom", fontproperties=CN_FONT, fontsize=7.2)
    fig.subplots_adjust(top=0.81, bottom=0.11, left=0.12, right=0.985, hspace=0.32)
    save_figure(fig, OUTPUT / "Fig_FE_path_local_50pct_S0_S90_cn")


def write_notes(metrics: pd.DataFrame, time: np.ndarray) -> None:
    dt_ms = 1000.0 * float(np.median(np.diff(time)))
    lines = [
        "0°/90°故障增量自身50%半峰值诊断",
        "===================================",
        "输入为保存的健康、0.90 mm太阳轮裂纹四测点FE路径差值与既有48个源事件；未重算动力学或FE模态。",
        "每次事件窗口为[源开始-4 ms,源结束+4 ms]，截断在记录范围内。",
        "峰值P_j=max|a_fault,j-a_healthy,j|（该事件窗口）；50%门槛T_j=0.5P_j。",
        "t50,j为从源事件起点起第一个满足|故障增量|>=T_j的采样时刻。",
        "两个测点分别使用自己的P_j，门槛绝对值不同；不能据此断言有限波速波前首达。",
        f"采样时间间隔={dt_ms:.9f} ms。",
        "48次事件由一套确定性仿真的四种重复几何相位组成，不是48次独立试验。",
    ]
    for state in (1, 2, 3, 4):
        group = metrics.loc[metrics["state"].eq(state)]
        n0 = int(group["earlier_half_peak_sensor_deg"].eq(0).sum())
        n90 = int(group["earlier_half_peak_sensor_deg"].eq(90).sum())
        mdt = float(group["first50_t90_minus_t0_ms"].median())
        n0amp = int(group["larger_peak_sensor_deg"].eq(0).sum())
        lines.append(
            f"状态{state}：半峰值时间0°先{n0}/{len(group)}、90°先{n90}/{len(group)}；"
            f"t90-t0中位数={mdt:+.6f} ms；0°幅值更大{n0amp}/{len(group)}。"
        )
    for number in (4, 8):
        row = metrics.loc[metrics["event_number"].eq(number)].iloc[0]
        lines.append(
            f"事件{number}：t50,0={1000*float(row['first50_sensor_0_lag_s']):.6f} ms；"
            f"t50,90={1000*float(row['first50_sensor_90_lag_s']):.6f} ms；"
            f"t90-t0={float(row['first50_t90_minus_t0_ms']):+.6f} ms。"
        )
    (OUTPUT / "threshold50_notes_cn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    events, time, values = load_inputs()
    metrics = compute_half_peak_metrics(events, time, values)
    metrics.to_csv(OUTPUT / "event_50pct_time_amplitude_0_90.csv", index=False,
                   encoding="utf-8-sig", quoting=csv.QUOTE_MINIMAL)
    state_figure(metrics)
    local_figure(metrics, time, values)
    write_notes(metrics, time)
    print(f"50% diagnostic saved: {OUTPUT}")
    print("earlier half-peak: 0°=", int(metrics["earlier_half_peak_sensor_deg"].eq(0).sum()),
          "90°=", int(metrics["earlier_half_peak_sensor_deg"].eq(90).sum()))


if __name__ == "__main__":
    main()
