#!/usr/bin/env python3
"""Build a 0°/90° evidence chain for the planet-crack FE transfer path.

This is a read-only replay of the saved 18-DOF healthy and planet-q=0.90 mm
records.  It uses the reflected PR-normal diagnostic asset, exact-ZOH modal
filters, and the same event windows as the exploratory path audit.  No source
force, crack parameter, modal shape, or gain is modified.

The output has three purposes:
1. quantify the pairwise distance trend (0° versus 90° only);
2. show one event that follows the distance trend and one boundary event;
3. demonstrate the modal explanation with full-window leave-one-mode-out
   diagnostics, rather than relying on a single instantaneous modal value.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from scipy.io import loadmat
from scipy.linalg import expm
from scipy.signal import lfilter, lfilter_zi
from scipy.stats import binomtest, pearsonr, spearmanr


PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "scripts"))
from figure_qa import require_matplotlib_panel_alignment  # noqa: E402


SOURCE = PROJECT / "results" / "thesis_ch2_fault_dynamics_20260912_094050"
ASSET = (
    PROJECT / "results" / "fe_planet_q090_23mode_exploratory_20260929"
    / "fe_asset_reflected_pr_normal_sensitivity.mat"
)
OUTPUT = PROJECT / "results" / "fe_planet_q090_23mode_exploratory_20260929"
SETTLED_START_S = 0.2
EVENT_TAIL_S = 0.004
ZETA = 0.02
PAIR_SENSOR_DEG = (0.0, 90.0)
PAIR_COLORS = {0.0: "#1f77b4", 90.0: "#c44e52"}


def read_state(path: Path) -> dict[str, np.ndarray | float]:
    raw = loadmat(path)
    fields = {
        "time": "time_s",
        "phase": "mesh_phase_cycles",
        "force": "mesh_force_pr_signed_N",
        "sp_loss": "fault_stiffness_loss_sp_N_per_m",
        "pr_loss": "fault_stiffness_loss_pr_N_per_m",
    }
    state: dict[str, np.ndarray | float] = {}
    for key, matlab_key in fields.items():
        if matlab_key not in raw:
            raise KeyError(f"{path} lacks {matlab_key}")
        state[key] = (
            np.asarray(raw[matlab_key], dtype=float).reshape(-1)
            if key in ("time", "phase")
            else np.asarray(raw[matlab_key], dtype=float)
        )
    state["fs"] = float(raw["fs_Hz"].item())
    state["dof"] = int(raw["model_DOF"].item())
    return state


def event_segments(time: np.ndarray, loss: np.ndarray, kind: str) -> list[dict[str, object]]:
    active = loss[:, 0] > max(float(np.max(loss[:, 0])) * 1e-10, 1e-6)
    edges = np.diff(np.r_[False, active, False].astype(np.int8))
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1) - 1
    result: list[dict[str, object]] = []
    for first, last in zip(starts, ends):
        if time[first] < SETTLED_START_S or last >= len(time) - 1:
            continue
        if time[last] + EVENT_TAIL_S > time[-1]:
            continue
        peak = first + int(np.argmax(loss[first:last + 1, 0]))
        window_end = int(np.searchsorted(time, time[last] + EVENT_TAIL_S, side="right"))
        result.append({"kind": kind, "first": int(first), "last": int(last), "peak": int(peak), "window": slice(int(first), window_end)})
    return result


def load_asset(path: Path) -> dict[str, np.ndarray]:
    raw = loadmat(path)
    asset = {
        "frequency_hz": np.asarray(raw["frequency_Hz"], dtype=float).reshape(-1),
        "contact_angle_rad": np.asarray(raw["contact_angle_rad"], dtype=float).reshape(-1),
        "input_shape_normal": np.asarray(raw["input_shape_normal"], dtype=float),
        "sensor_angle_local_deg": np.asarray(raw["sensor_angle_local_deg"], dtype=float).reshape(-1),
        "sensor_angle_global_deg": np.asarray(raw["sensor_angle_global_deg"], dtype=float).reshape(-1),
        "sensor_shape_radial": np.asarray(raw["sensor_shape_radial"], dtype=float),
    }
    order = np.argsort(asset["contact_angle_rad"])
    asset["contact_angle_rad"] = asset["contact_angle_rad"][order]
    asset["input_shape_normal"] = asset["input_shape_normal"][order, :]
    return asset


def contact_angles(carrier_angle_rad: np.ndarray, sensor_global_deg: np.ndarray) -> np.ndarray:
    phases = np.deg2rad(np.array([0.0, 120.0, 240.0]))[None, :]
    carrier = np.asarray(carrier_angle_rad, dtype=float)[:, None]
    return np.mod(np.deg2rad(float(sensor_global_deg[0])) - carrier - phases, 2.0 * np.pi)


def exact_zoh_acceleration_filter(omega: float, zeta: float, dt: float) -> tuple[np.ndarray, np.ndarray]:
    a_cont = np.array([[0.0, 1.0], [-omega * omega, -2.0 * zeta * omega]])
    augmented = np.zeros((3, 3), dtype=float)
    augmented[:2, :2] = a_cont
    augmented[:2, 2] = (0.0, 1.0)
    matrix = expm(augmented * dt)
    ad, bd = matrix[:2, :2], matrix[:2, 2]
    c = np.array([-omega * omega, -2.0 * zeta * omega])
    trace, det = float(np.trace(ad)), float(np.linalg.det(ad))
    a = np.array([1.0, -trace, det])
    b = np.array([1.0, -trace + float(c @ bd), det + float(c @ (ad - trace * np.eye(2)) @ bd)])
    return b, a


def save_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build_components(
    time: np.ndarray,
    delta_force: np.ndarray,
    carrier: np.ndarray,
    asset: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    """Return n×2×3×nm signed acceleration components and FE angles."""

    fs = 1.0 / float(np.median(np.diff(time)))
    frequency = np.asarray(asset["frequency_hz"], dtype=float)
    theta = contact_angles(carrier, asset["sensor_angle_global_deg"])
    angle_grid = np.asarray(asset["contact_angle_rad"], dtype=float)
    shape_in = np.asarray(asset["input_shape_normal"], dtype=float)
    sensor_local = np.asarray(asset["sensor_angle_local_deg"], dtype=float)
    pair_idx = [int(np.flatnonzero(np.isclose(sensor_local, x))[0]) for x in PAIR_SENSOR_DEG]
    shape_sensor = np.asarray(asset["sensor_shape_radial"], dtype=float)[pair_idx, :]
    n, nm = time.size, frequency.size
    component = np.empty((n, 2, 3, nm), dtype=np.float64)
    dt = 1.0 / fs
    xp = np.r_[angle_grid[-1] - 2.0 * np.pi, angle_grid, angle_grid[0] + 2.0 * np.pi]
    if not np.all(np.diff(xp) > 0.0):
        raise ValueError("Periodic contact-angle interpolation grid must be strictly increasing.")
    for im, hz in enumerate(frequency):
        fp = np.r_[shape_in[-1, im], shape_in[:, im], shape_in[0, im]]
        b, a = exact_zoh_acceleration_filter(2.0 * np.pi * hz, ZETA, dt)
        for ip in range(3):
            psi = np.interp(theta[:, ip], xp, fp)
            q = delta_force[:, ip] * psi
            eta_ddot, _ = lfilter(b, a, q, zi=lfilter_zi(b, a) * q[0])
            component[:, :, ip, im] = eta_ddot[:, None] * shape_sensor[:, im][None, :]
    return component, theta


def make_event_summary(
    time: np.ndarray,
    component: np.ndarray,
    theta: np.ndarray,
    events: list[dict[str, object]],
    sensor_global: np.ndarray,
    delta_force: np.ndarray,
) -> tuple[list[dict[str, object]], dict[int, dict[str, object]]]:
    total = component.sum(axis=(2, 3))
    rows: list[dict[str, object]] = []
    event_data: dict[int, dict[str, object]] = {}
    for number, event in enumerate(events, start=1):
        window = event["window"]
        assert isinstance(window, slice)
        first, last, peak = (int(event[key]) for key in ("first", "last", "peak"))
        contact_global = float(np.rad2deg(theta[peak, 0]) % 360.0)
        contact_local = float((sensor_global[0] - contact_global) % 360.0)
        distances = np.abs((contact_local - np.asarray(PAIR_SENSOR_DEG) + 180.0) % 360.0 - 180.0)
        p0 = float(np.max(np.abs(total[window, 0])))
        p90 = float(np.max(np.abs(total[window, 1])))
        peaks = np.array([p0, p90])
        near_idx = int(np.argmin(distances))
        far_idx = 1 - near_idx
        near_peak, far_peak = float(peaks[near_idx]), float(peaks[far_idx])
        near_label, far_label = PAIR_SENSOR_DEG[near_idx], PAIR_SENSOR_DEG[far_idx]
        first50 = []
        for isensor in range(2):
            values = np.abs(total[window, isensor])
            threshold = 0.5 * float(np.max(values))
            cross = np.flatnonzero(values >= threshold)
            first50.append(float(time[window.start + int(cross[0])]))
        force_peak = float(np.max(np.abs(delta_force[window, 0])))
        row = {
            "event_number": number,
            "source_mesh": str(event["kind"]),
            "contact_position_deg": contact_local,
            "distance_sensor_0_deg": float(distances[0]),
            "distance_sensor_90_deg": float(distances[1]),
            "near_sensor_deg": near_label,
            "far_sensor_deg": far_label,
            "peak_sensor_0_ms2": p0,
            "peak_sensor_90_ms2": p90,
            "peak_sensor_0_scaled_1e3": p0 * 1e3,
            "peak_sensor_90_scaled_1e3": p90 * 1e3,
            "peak_near_ms2": near_peak,
            "peak_far_ms2": far_peak,
            "near_far_peak_ratio": near_peak / far_peak,
            "winner_sensor_deg": PAIR_SENSOR_DEG[int(np.argmax(peaks))],
            "near_equals_winner": int(near_idx == int(np.argmax(peaks))),
            "first50_sensor_0_s": first50[0],
            "first50_sensor_90_s": first50[1],
            "near_minus_far_first50_us": (first50[near_idx] - first50[far_idx]) * 1e6,
            "peak_per_pr1_delta_0_scaled_1e3": p0 / force_peak * 1e3,
            "peak_per_pr1_delta_90_scaled_1e3": p90 / force_peak * 1e3,
            "source_peak_time_s": float(time[peak]),
            "source_start_s": float(time[first]),
            "source_end_s": float(time[last]),
        }
        rows.append(row)
        event_data[number] = {
            "event": event,
            "window": window,
            "contact_position_deg": contact_local,
            "total": total[window, :].copy(),
            "component": component[window, :, :, :].copy(),
            "time": time[window].copy(),
            "summary": row,
        }
    return rows, event_data


def modal_metrics(event_data: dict[int, dict[str, object]], frequencies: np.ndarray) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for event_number, data in event_data.items():
        component = np.asarray(data["component"])
        total = np.asarray(data["total"])
        time = np.asarray(data["time"])
        for isensor, sensor in enumerate(PAIR_SENSOR_DEG):
            mode_signals = component[:, isensor, :, :].sum(axis=1)  # time × mode
            net = mode_signals.sum(axis=1)
            peak_index = int(np.argmax(np.abs(net)))
            peak_abs = float(np.max(np.abs(net)))
            rms_total = float(np.sqrt(np.mean(net * net)))
            rms_modes = np.sqrt(np.mean(mode_signals * mode_signals, axis=0))
            abs_rms_sum = float(np.sum(np.abs(rms_modes)))
            for im, hz in enumerate(frequencies):
                mode = mode_signals[:, im]
                without = net - mode
                peak_without = float(np.max(np.abs(without)))
                rms_without = float(np.sqrt(np.mean(without * without)))
                fixed_peak_without = abs(float(net[peak_index] - mode[peak_index]))
                rows.append({
                    "event_number": event_number,
                    "sensor_local_deg": sensor,
                    "mode_index": im + 1,
                    "mode_frequency_hz": float(hz),
                    "mode_signed_at_total_peak_ms2": float(mode[peak_index]),
                    "mode_peak_abs_ms2": float(np.max(np.abs(mode))),
                    "mode_rms_ms2": float(rms_modes[im]),
                    "absolute_rms_share": float(abs(rms_modes[im]) / abs_rms_sum) if abs_rms_sum else 0.0,
                    "total_peak_abs_ms2": peak_abs,
                    "total_rms_ms2": rms_total,
                    "peak_without_mode_abs_ms2": peak_without,
                    "peak_drop_full_window_fraction": float((peak_abs - peak_without) / peak_abs),
                    "peak_drop_at_fixed_total_peak_fraction": float((peak_abs - fixed_peak_without) / peak_abs),
                    "rms_drop_full_window_fraction": float((rms_total - rms_without) / rms_total),
                    "total_peak_time_s": float(time[peak_index]),
                })
    return rows


def get_mode_rows(rows: list[dict[str, object]], event: int, sensor: float) -> list[dict[str, object]]:
    return [r for r in rows if int(r["event_number"]) == event and float(r["sensor_local_deg"]) == sensor]


def configure_fonts() -> tuple[FontProperties, FontProperties]:
    mpl.rcParams.update({
        "font.family": ["Times New Roman", "SimSun", "Arial"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 8.2,
        "axes.labelsize": 8.2,
        "axes.titlesize": 8.8,
        "axes.linewidth": 0.65,
        "xtick.labelsize": 7.4,
        "ytick.labelsize": 7.4,
        "legend.fontsize": 7.4,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })
    cn = FontProperties(family="SimSun", size=8.4)
    cn_small = FontProperties(family="SimSun", size=7.4)
    return cn, cn_small


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.14, 1.04, label, transform=ax.transAxes, fontweight="bold", fontsize=9, va="bottom")


def plot_evidence(summary_rows: list[dict[str, object]], modal_rows: list[dict[str, object]], event_data: dict[int, dict[str, object]], output: Path) -> None:
    cn, cn_small = configure_fonts()
    cn_unicode = FontProperties(family="Arial Unicode MS", size=8.4)
    fig, axes = plt.subplots(3, 2, figsize=(7.25, 9.4), constrained_layout=False)
    axa, axb, axc, axd, axe, axf = axes.ravel()
    # (a) paired spatial trend
    for row in summary_rows:
        x = float(row["contact_position_deg"])
        y0, y90 = float(row["peak_per_pr1_delta_0_scaled_1e3"]), float(row["peak_per_pr1_delta_90_scaled_1e3"])
        highlighted = int(row["event_number"]) in (7, 9)
        marker_size = 39 if highlighted else 26
        edge_width = 0.85 if highlighted else 0.3
        axa.plot([x, x], [y0, y90], color="#b8b8b8", linewidth=0.65, zorder=1)
        axa.scatter(x, y0, s=marker_size, color=PAIR_COLORS[0.0], edgecolor="black", linewidth=edge_width, zorder=3)
        axa.scatter(x, y90, s=marker_size, color=PAIR_COLORS[90.0], edgecolor="black", linewidth=edge_width, zorder=3)
    axa.axvline(0, color=PAIR_COLORS[0.0], alpha=0.45, linewidth=0.8)
    axa.axvline(90, color=PAIR_COLORS[90.0], alpha=0.45, linewidth=0.8)
    axa.set_xlim(-8, 368)
    axa.set_xlabel("PR 接触位置 / (°)", fontproperties=cn, labelpad=3)
    axa.set_ylabel("单位故障力峰值 / [10⁻³ (m/s²)/N]", fontproperties=cn_unicode, labelpad=3)
    axa.set_title("0°/90°测点：响应随啮合位置变化", fontproperties=cn, pad=5)
    axa.grid(axis="y", color="#dedede", linewidth=0.45)
    # (b) near/far paired comparison
    events = [int(r["event_number"]) for r in summary_rows]
    ratios = [float(r["near_far_peak_ratio"]) for r in summary_rows]
    colors = ["#2a9d8f" if x >= 1 else "#d9822b" for x in ratios]
    axb.bar(events, ratios, color=colors, edgecolor="#333333", linewidth=0.35, width=0.68)
    axb.axhline(1.0, color="#222222", linewidth=0.8)
    axb.set_xticks(events)
    axb.set_xlabel("事件编号", fontproperties=cn, labelpad=3)
    axb.set_ylabel("近/远测点峰值比", fontproperties=cn, labelpad=3)
    axb.set_title("10次事件中，近测点7次峰值更大", fontproperties=cn, pad=5)
    axb.grid(axis="y", color="#dedede", linewidth=0.45)
    # waveform helper
    def waveform(ax: plt.Axes, event_number: int, title: str) -> None:
        data = event_data[event_number]
        t = (np.asarray(data["time"]) - float(data["time"][0])) * 1000.0
        total = np.asarray(data["total"]) * 1e3
        for i, sensor in enumerate(PAIR_SENSOR_DEG):
            ax.plot(t, total[:, i], color=PAIR_COLORS[sensor], linewidth=0.95, label=f"Sensor {int(sensor)}°")
            idx = int(np.argmax(np.abs(total[:, i])))
            ax.scatter(t[idx], total[idx, i], color=PAIR_COLORS[sensor], edgecolor="black", linewidth=0.3, s=22, zorder=4)
        ax.axhline(0, color="#555555", linewidth=0.45)
        ax.set_xlabel("事件起点后的时间 / ms", fontproperties=cn, labelpad=3)
        ax.set_ylabel("故障增量 / (10⁻³ m/s²)", fontproperties=cn_unicode, labelpad=3)
        ax.set_title(title, fontproperties=cn, pad=5)
        ax.grid(axis="y", color="#dedede", linewidth=0.45)
    waveform(axc, 9, "事件9：近测点0°峰值更大")
    waveform(axe, 7, "事件7：近测点0°但90°峰值更大")
    # modal leave-one-out helper
    def modal_panel(ax: plt.Axes, event_number: int, title: str) -> None:
        selected = []
        for sensor in PAIR_SENSOR_DEG:
            rows = get_mode_rows(modal_rows, event_number, sensor)
            selected.extend(rows)
        by_mode = {}
        for row in selected:
            k = int(row["mode_index"])
            by_mode[k] = max(by_mode.get(k, -np.inf), float(row["peak_drop_full_window_fraction"]))
        top_modes = sorted(by_mode, key=by_mode.get, reverse=True)[:6]
        top_modes = list(reversed(top_modes))
        y = np.arange(len(top_modes))
        height = 0.34
        for offset, sensor in [(-height / 2, 0.0), (height / 2, 90.0)]:
            vals = []
            labels = []
            for k in top_modes:
                rr = [r for r in modal_rows if int(r["event_number"]) == event_number and float(r["sensor_local_deg"]) == sensor and int(r["mode_index"]) == k][0]
                vals.append(float(rr["peak_drop_full_window_fraction"]) * 100.0)
                labels.append(f"{k} ({float(rr['mode_frequency_hz'])/1000:.2f})")
            ax.barh(y + offset, vals, height=height, color=PAIR_COLORS[sensor], edgecolor="#333333", linewidth=0.3, label=f"Sensor {int(sensor)}°")
        ax.set_yticks(y, labels)
        ax.set_xlabel("剔除单阶模态后的峰值下降 / %", fontproperties=cn, labelpad=3)
        ax.set_ylabel("模态阶次（频率 / kHz）", fontproperties=cn, labelpad=3)
        ax.axvline(0, color="#333333", linewidth=0.65)
        ax.set_title(title, fontproperties=cn, pad=5)
        ax.grid(axis="x", color="#dedede", linewidth=0.45)
    modal_panel(axd, 9, "事件9：逐阶剔除后的峰值敏感性")
    modal_panel(axf, 7, "事件7：逐阶剔除后的峰值敏感性")
    for ax, label in zip((axa, axb, axc, axd, axe, axf), "abcdef"):
        add_panel_label(ax, f"({label})")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("行星轮裂纹传递路径证据链：0°与90°测点", fontproperties=cn, fontsize=11, y=0.985)
    fig.text(0.5, 0.012, "蓝/红：0°/90°测点；粗边圆点：事件7/9；绿/橙：近测点占优/未占优；响应为故障减健康；模态面板为整段事件窗逐阶剔除结果。", ha="center", fontproperties=cn_small, color="#444444")
    fig.subplots_adjust(left=0.105, right=0.975, bottom=0.075, top=0.935, hspace=0.62, wspace=0.36)
    require_matplotlib_panel_alignment(
        fig,
        json_out=str(output) + ".alignment.json",
        overlay_svg=str(output) + ".alignment.svg",
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        strict=True,
    )
    fig.savefig(str(output) + ".png", dpi=600, bbox_inches="tight")
    fig.savefig(str(output) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(output) + ".pdf", bbox_inches="tight")
    fig.savefig(str(output) + ".svg", bbox_inches="tight")
    plt.close(fig)


def plot_geometry_timing(summary_rows: list[dict[str, object]], output: Path) -> None:
    """Show the local angle convention, event order, and threshold timing."""

    cn, cn_small = configure_fonts()
    fig = plt.figure(figsize=(7.25, 3.25), constrained_layout=False)
    grid = fig.add_gridspec(1, 3, wspace=0.42)
    axa = fig.add_subplot(grid[0, 0], projection="polar")
    axb = fig.add_subplot(grid[0, 1])
    axc = fig.add_subplot(grid[0, 2])
    events = np.array([int(row["event_number"]) for row in summary_rows])
    angles_deg = np.array([float(row["contact_position_deg"]) for row in summary_rows])
    source_mesh = [str(row["source_mesh"]) for row in summary_rows]
    lead_us = -np.array([float(row["near_minus_far_first50_us"]) for row in summary_rows])
    mesh_color = {"PR": "#2a9d8f", "SP": "#d9822b"}
    mesh_marker = {"PR": "o", "SP": "s"}

    # (a) Local angular convention: 0° at the top and clockwise positive.
    axa.set_theta_zero_location("N")
    axa.set_theta_direction(-1)
    axa.set_ylim(0.0, 1.14)
    theta = np.deg2rad(angles_deg)
    for number, angle, mesh in zip(events, theta, source_mesh):
        axa.scatter(angle, 0.82, s=28, marker=mesh_marker[mesh], color=mesh_color[mesh],
                    edgecolor="black", linewidth=0.35, zorder=3)
        axa.text(angle, 0.63, str(number), ha="center", va="center", fontsize=7.2)
    axa.scatter(0.0, 1.02, marker="^", s=52, color=PAIR_COLORS[0.0], edgecolor="black", linewidth=0.45, zorder=4)
    axa.scatter(np.deg2rad(90.0), 1.02, marker="^", s=52, color=PAIR_COLORS[90.0], edgecolor="black", linewidth=0.45, zorder=4)
    axa.set_xticks(np.deg2rad([0, 90, 180, 270]), ["0°", "90°", "180°", "270°"])
    axa.set_yticks([])
    axa.grid(False)
    axa.spines["polar"].set_color("#555555")
    axa.spines["polar"].set_linewidth(0.65)
    axa.set_title("本地角度与10次事件位置", fontproperties=cn, pad=8)

    # (b) Event order and the concurrent P1 PR input position on the ring.
    axb.plot(events, angles_deg, color="#888888", linewidth=0.7, zorder=1)
    for number, angle, mesh in zip(events, angles_deg, source_mesh):
        axb.scatter(number, angle, s=31, marker=mesh_marker[mesh], color=mesh_color[mesh],
                    edgecolor="black", linewidth=0.35, zorder=3)
    axb.axhline(0.0, color=PAIR_COLORS[0.0], linewidth=0.8, alpha=0.75)
    axb.axhline(90.0, color=PAIR_COLORS[90.0], linewidth=0.8, alpha=0.75)
    axb.set_xlim(0.5, 10.5)
    axb.set_ylim(-12.0, 372.0)
    axb.set_xticks(events)
    axb.set_yticks([0, 90, 180, 270, 360])
    axb.set_xlabel("事件编号", fontproperties=cn, labelpad=3)
    axb.set_ylabel("P1 的 PR 接触位置 / (°)", fontproperties=cn, labelpad=3)
    axb.set_title("事件序列中的啮合位置", fontproperties=cn, pad=8)
    axb.grid(axis="y", color="#dedede", linewidth=0.45)

    # (c) Timing metric.  Positive means the nearer member of the 0°/90° pair
    # reaches 50% of its own event peak first.
    axc.bar(events, lead_us, color="#4c78a8", edgecolor="#333333", linewidth=0.35, width=0.68)
    axc.axhline(0.0, color="#333333", linewidth=0.65)
    axc.set_xlim(0.5, 10.5)
    axc.set_xticks(events)
    axc.set_xlabel("事件编号", fontproperties=cn, labelpad=3)
    axc.set_ylabel("近测点50%阈值领先量 / μs", fontproperties=cn, labelpad=3)
    axc.set_title("10/10次近测点先越过自身50%阈值", fontproperties=cn, pad=8)
    axc.grid(axis="y", color="#dedede", linewidth=0.45)

    for ax, label in zip((axa, axb, axc), "abc"):
        add_panel_label(ax, f"({label})")
        if ax is not axa:
            ax.set_box_aspect(1.0)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
    fig.suptitle("行星轮裂纹事件的啮合位置与响应先后", fontproperties=cn, fontsize=11, y=0.975)
    fig.text(
        0.5,
        0.018,
        "蓝/红三角：0°/90°测点；圆/方：PR/SP故障事件；角度均为齿圈上的P1 PR输入位置。50%阈值仅表示波形先后，不等于材料波传播时间。",
        ha="center",
        fontproperties=cn_small,
        color="#444444",
    )
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.20, top=0.82, wspace=0.46)
    require_matplotlib_panel_alignment(
        fig,
        json_out=str(output) + ".alignment.json",
        overlay_svg=str(output) + ".alignment.svg",
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        strict=True,
    )
    fig.savefig(str(output) + ".png", dpi=600, bbox_inches="tight")
    fig.savefig(str(output) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(output) + ".pdf", bbox_inches="tight")
    fig.savefig(str(output) + ".svg", bbox_inches="tight")
    plt.close(fig)


def write_report(path: Path, summary_rows: list[dict[str, object]], modal_rows: list[dict[str, object]], event_data: dict[int, dict[str, object]], reconstruction_error: float) -> dict[str, object]:
    # Use all 20 paired observations only as a descriptive distance diagnostic.
    # The two sensor values within each event are paired and the ten events come
    # from one deterministic simulation, so p-values are not treated as
    # inferential evidence.
    d_all, p_norm_all, p_raw_all = [], [], []
    for row in summary_rows:
        d_all.extend([float(row["distance_sensor_0_deg"]), float(row["distance_sensor_90_deg"])])
        p_norm_all.extend([float(row["peak_per_pr1_delta_0_scaled_1e3"]), float(row["peak_per_pr1_delta_90_scaled_1e3"])])
        p_raw_all.extend([float(row["peak_sensor_0_ms2"]), float(row["peak_sensor_90_ms2"])])
    pearson = pearsonr(d_all, p_norm_all)
    spearman = spearmanr(d_all, p_norm_all)
    pearson_raw = pearsonr(d_all, p_raw_all)
    spearman_raw = spearmanr(d_all, p_raw_all)
    near_wins = sum(int(r["near_equals_winner"]) for r in summary_rows)
    earlier = sum(float(r["near_minus_far_first50_us"]) < 0 for r in summary_rows)
    pair_sign_test = binomtest(near_wins, len(summary_rows), p=0.5, alternative="two-sided")
    signed_first50 = [float(r["near_minus_far_first50_us"]) for r in summary_rows]
    positive_lead = [-value for value in signed_first50]
    report_rows = []
    event_details = []
    for event_number in (9, 7):
        summary = next(row for row in summary_rows if int(row["event_number"]) == event_number)
        cancellation = {}
        for sensor in PAIR_SENSOR_DEG:
            rr = get_mode_rows(modal_rows, event_number, sensor)
            top = sorted(rr, key=lambda x: float(x["peak_drop_full_window_fraction"]), reverse=True)[:5]
            total_peak = float(rr[0]["total_peak_abs_ms2"])
            cancellation[str(int(sensor))] = float(
                sum(abs(float(x["mode_signed_at_total_peak_ms2"])) for x in rr) / total_peak
            )
            report_rows.append({
                "event": event_number,
                "sensor": sensor,
                "top_modes_by_full_window_peak_drop": [(int(x["mode_index"]), round(float(x["mode_frequency_hz"]), 2), round(float(x["peak_drop_full_window_fraction"])*100, 2)) for x in top],
            })
        event_details.append({
            "event": event_number,
            "source_mesh": str(summary["source_mesh"]),
            "contact_position_deg": float(summary["contact_position_deg"]),
            "distance_sensor_0_deg": float(summary["distance_sensor_0_deg"]),
            "distance_sensor_90_deg": float(summary["distance_sensor_90_deg"]),
            "peak_sensor_0_ms2": float(summary["peak_sensor_0_ms2"]),
            "peak_sensor_90_ms2": float(summary["peak_sensor_90_ms2"]),
            "peak_0_over_90": float(summary["peak_sensor_0_ms2"]) / float(summary["peak_sensor_90_ms2"]),
            "pair_near_sensor_deg": float(summary["near_sensor_deg"]),
            "pair_winner_sensor_deg": float(summary["winner_sensor_deg"]),
            "cancellation_coefficient_at_total_peak": cancellation,
        })
    result = {
        "scope": "0° and 90° sensors only; 120° and 240° excluded from the primary comparison",
        "claim": "Angular proximity produces a repeatable tendency rather than a strict rule; modal projection and signed superposition explain event-specific deviations.",
        "events": 10,
        "near_sensor_peak_wins": near_wins,
        "near_sensor_peak_wins_sign_test_p_descriptive_only": float(pair_sign_test.pvalue),
        "near_sensor_first_50pct_earlier": earlier,
        "near_far_ratio_median": float(np.median([float(r["near_far_peak_ratio"]) for r in summary_rows])),
        "near_far_ratio_mean": float(np.mean([float(r["near_far_peak_ratio"]) for r in summary_rows])),
        "distance_force_normalized_peak_pearson_r": float(pearson.statistic),
        "distance_force_normalized_peak_pearson_p_descriptive_only": float(pearson.pvalue),
        "distance_force_normalized_peak_spearman_rho": float(spearman.statistic),
        "distance_force_normalized_peak_spearman_p_descriptive_only": float(spearman.pvalue),
        "distance_raw_peak_pearson_r": float(pearson_raw.statistic),
        "distance_raw_peak_spearman_rho": float(spearman_raw.statistic),
        "near_minus_far_first50_us_range": [float(min(signed_first50)), float(max(signed_first50))],
        "near_first50_lead_us_range": [float(min(positive_lead)), float(max(positive_lead))],
        "component_reconstruction_relative_l2": reconstruction_error,
        "representative_events": report_rows,
        "representative_event_details": event_details,
        "statistical_boundary": "The correlations describe 20 paired sensor observations from ten events in one deterministic run; their p-values are not used as independent-sample significance evidence.",
        "interpretation_boundary": "The modal evidence is conditional on the reflected-normal 23-mode diagnostic replay; it does not validate SI normalization, omitted-mode convergence, or a travelling-wave speed.",
    }
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    e9 = next(item for item in event_details if int(item["event"]) == 9)
    e7 = next(item for item in event_details if int(item["event"]) == 7)
    e9_mode9 = next(
        row for row in get_mode_rows(modal_rows, 9, 0.0)
        if int(row["mode_index"]) == 9
    )
    e7_mode1_90 = next(
        row for row in get_mode_rows(modal_rows, 7, 90.0)
        if int(row["mode_index"]) == 1
    )
    txt = [
        "行星轮裂纹齿圈传递路径证据链（0°/90°主比较）",
        "",
        "结论：角距离呈可重复的总体趋势，但不是严格单调定律；齿圈模态输入投影、测点振型可观测性及有符号模态叠加会造成事件级偏离。",
        f"完整事件：{len(summary_rows)}次；近测点峰值更大：{near_wins}/{len(summary_rows)}；近测点50%阈值更早：{earlier}/{len(summary_rows)}。",
        f"7/10仅作为模型内描述；即使把10次事件视为独立，双侧符号检验p={result['near_sensor_peak_wins_sign_test_p_descriptive_only']:.3f}，不足以声称严格距离定律。",
        f"近/远峰值比：均值{result['near_far_ratio_mean']:.2f}，中位数{result['near_far_ratio_median']:.2f}。",
        f"与图(a)一致的20个力归一化观测中，距离—峰值 Pearson r={result['distance_force_normalized_peak_pearson_r']:.3f}，Spearman rho={result['distance_force_normalized_peak_spearman_rho']:.3f}；仅作描述，不作独立样本显著性推断。",
        f"近测点越过自身50%阈值的领先量为{result['near_first50_lead_us_range'][0]:.1f}–{result['near_first50_lead_us_range'][1]:.1f} μs；该指标是波形阈值，不是传播速度。",
        "",
        "证据1：事件位置—0°/90°峰值配对图显示总体随角距离增大而衰减。",
        "证据2：近/远峰值比图给出7/10次近测点占优，并显式保留3个偏离事件。",
        "证据3：事件9是距离规律成立的例子（近0°且0°更强）；事件7是边界例子（近0°但90°更强）。",
        "证据4：每阶模态响应在整段事件窗口内逐阶剔除，比较去掉该模态后的峰值下降，而不是只看一个采样时刻。",
        f"模态分项相加重构误差：{reconstruction_error:.3e}，说明分解在数值上闭合。",
        "",
        f"事件9：接触位置{e9['contact_position_deg']:.3f}°，距0°/90°测点分别{e9['distance_sensor_0_deg']:.3f}°/{e9['distance_sensor_90_deg']:.3f}°；0°峰值为90°的{e9['peak_0_over_90']:.3f}倍。第9阶（19.254 kHz）剔除后，0°窗口峰值下降{float(e9_mode9['peak_drop_full_window_fraction'])*100:.2f}%，RMS下降{float(e9_mode9['rms_drop_full_window_fraction'])*100:.2f}%。",
        f"事件7：接触位置{e7['contact_position_deg']:.3f}°；在0°/90°二测点内0°更近，但90°峰值为0°的{1.0/float(e7['peak_0_over_90']):.3f}倍。总峰值时刻的相消系数为0°={e7['cancellation_coefficient_at_total_peak']['0']:.3f}、90°={e7['cancellation_coefficient_at_total_peak']['90']:.3f}；90°第1阶（15.104 kHz）瞬时贡献为{float(e7_mode1_90['mode_signed_at_total_peak_ms2']):.6f} m/s²。",
        "逐阶剔除指标若为负，表示去掉该模态后整段窗口的最大峰值反而增大，反映相消关系或峰值时刻迁移；不能把该百分数解释成方差占比。",
        "",
        "模态解释：Q_r(t)=Σ_i ΔF_PR,i(t)·ψ_in,r[θ_i(t)]；a_s(t)=Σ_r φ_out,s,r·η̈_r(t)。啮合位置决定输入振型投影，测点位置决定输出振型可观测性，模态频率/阻尼决定动态放大，最终各模态有符号叠加。",
        "角距离并未作为经验增益或人工延时写入方程；距离规律是空间振型共同作用后呈现的统计倾向。",
        "因此距离只提供总体倾向，不足以单独预测每次峰值；事件级偏离是模态投影和相位叠加的结果。",
        "",
        "数据版本：本图只使用反射PR法向的23阶诊断回放结果，不与同目录旧的原始法向统计混用。",
        "边界：本证据链使用51.2 kHz源采样可保留的23/60阶模态；ANSYS SI质量归一化、PR载荷方向和全模态截断收敛仍需独立核验。",
    ]
    (path.parent / "planet_modal_evidence_chain_cn.txt").write_text("\n".join(txt) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--asset", type=Path, default=ASSET)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = parser.parse_args()
    healthy = read_state(args.source / "case_healthy_state.mat")
    fault = read_state(args.source / "case_planet_root_crack_qc0p900mm_state.mat")
    time = np.asarray(fault["time"])
    # The saved source mesh counter is Zr*carrier_angle/(2*pi), with Zr=84.
    carrier = 2.0 * np.pi * np.asarray(fault["phase"]) / 84.0
    if not np.array_equal(time, np.asarray(healthy["time"])):
        raise ValueError("Healthy and fault time grids differ.")
    if int(fault["dof"]) != 18 or int(healthy["dof"]) != 18:
        raise ValueError("Expected 18-DOF source records.")
    delta_force = np.asarray(fault["force"]) - np.asarray(healthy["force"])
    asset = load_asset(args.asset)
    all_frequency = np.asarray(asset["frequency_hz"])
    keep = all_frequency < float(fault["fs"]) / 2.0
    asset["frequency_hz"] = all_frequency[keep]
    asset["input_shape_normal"] = np.asarray(asset["input_shape_normal"])[:, keep]
    asset["sensor_shape_radial"] = np.asarray(asset["sensor_shape_radial"])[:, keep]
    theta = contact_angles(carrier, asset["sensor_angle_global_deg"])
    sp_loss = np.asarray(fault["sp_loss"])
    pr_loss = np.asarray(fault["pr_loss"])
    events = event_segments(time, sp_loss, "SP") + event_segments(time, pr_loss, "PR")
    events.sort(key=lambda event: time[int(event["peak"])])
    if len(events) != 10:
        raise ValueError(f"Expected 10 complete events, got {len(events)}")
    component, theta = build_components(time, delta_force, carrier, asset)
    summary_rows, event_data = make_event_summary(time, component, theta, events, asset["sensor_angle_global_deg"], delta_force)
    modal_rows = modal_metrics(event_data, asset["frequency_hz"])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    save_csv(args.output_dir / "planet_pair_event_summary_0_90.csv", summary_rows)
    save_csv(args.output_dir / "planet_modal_window_metrics_0_90.csv", modal_rows)
    total = component.sum(axis=(2, 3))
    replay = np.zeros_like(total)
    for im, hz in enumerate(asset["frequency_hz"]):
        mode_total = component[:, :, :, im].sum(axis=2)
        replay += mode_total
    reconstruction_error = float(np.linalg.norm(replay - total) / np.linalg.norm(replay))
    plot_stem = args.output_dir / "Fig_planet_q090_modal_evidence_0_90"
    plot_evidence(summary_rows, modal_rows, event_data, plot_stem)
    plot_geometry_timing(
        summary_rows,
        args.output_dir / "Fig_planet_q090_geometry_timing_0_90_cn",
    )
    report = write_report(args.output_dir / "planet_modal_evidence_chain.json", summary_rows, modal_rows, event_data, reconstruction_error)
    print(json.dumps({k: report[k] for k in ("near_sensor_peak_wins", "near_sensor_first_50pct_earlier", "near_far_ratio_median", "distance_force_normalized_peak_pearson_r", "distance_force_normalized_peak_spearman_rho", "component_reconstruction_relative_l2")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
