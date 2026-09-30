#!/usr/bin/env python3
"""Plot evidence for the registered FE path and event-wise phase ordering.

The figures are descriptive post-processing of the saved 18-DOF source record
and the corrected four-sensor FE modal projection.  No waveform is shifted,
scaled, smoothed, or regenerated.  In particular, the phase figure uses the
original event clock and marks modal-ringdown peak times; it does not claim a
finite-wave-speed arrival time.
"""

from __future__ import annotations

import csv
import sys
from collections import OrderedDict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))
from figure_qa import require_matplotlib_panel_alignment


RESULTS = ROOT / "results" / "fe_causal_q090_corrected_registration_20260928"
EVENT_CSV = RESULTS / "corrected_event_metrics.csv"
SIGNAL_CSV = RESULTS / "corrected_registration_time_signals.csv"
OUT = RESULTS / "mapping_evidence"

SENSOR_ORDER = [0.0, 90.0, 120.0, 240.0]
SENSOR_COLORS = {
    0.0: "#1769aa",
    90.0: "#d95f02",
    120.0: "#2a9d55",
    240.0: "#7b4ab5",
}
STATE_ORDER = OrderedDict(
    [
        (4, ("A", 87.2241209)),
        (1, ("B", 177.2241209)),
        (2, ("C", 267.2259520)),
        (3, ("D", 357.2222899)),
    ]
)

mpl.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 8.2,
        "axes.titlesize": 9,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.2,
        "ytick.labelsize": 7.2,
        "legend.fontsize": 7.3,
        "axes.linewidth": 0.75,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    }
)


def save_figure(fig: mpl.figure.Figure, stem: Path, multi_panel: bool) -> None:
    """Run the required layout gate, then export editable/vector and raster files."""

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


def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    events = pd.read_csv(EVENT_CSV)
    signals = pd.read_csv(SIGNAL_CSV)
    events["event_number"] = events["event_number"].astype(int)
    events["state"] = (events["event_number"] - 1) % 4 + 1
    for col in (
        "sensor_angle_local_deg",
        "sensor_angle_global_deg",
        "contact_angle_global_deg",
        "contact_sensor_angular_distance_deg",
        "increment_peak_abs_ms2",
        "peak_lag_s",
        "first_5pct_peak_crossing_lag_s",
    ):
        events[col] = events[col].astype(float)
    return events, signals


def event_summary(events: pd.DataFrame) -> pd.DataFrame:
    """Make one row per source event, preserving all measured definitions."""

    rows: list[dict[str, object]] = []
    for event_number, group in events.groupby("event_number", sort=True):
        group = group.sort_values("sensor_angle_local_deg")
        first = group.iloc[0]
        nearest = group.loc[group["contact_sensor_angular_distance_deg"].idxmin()]
        strongest = group.loc[group["increment_peak_abs_ms2"].idxmax()]
        s0 = group.loc[group["sensor_angle_local_deg"].eq(0.0)].iloc[0]
        s90 = group.loc[group["sensor_angle_local_deg"].eq(90.0)].iloc[0]
        near_0_90 = 0.0 if float(s0["contact_sensor_angular_distance_deg"]) < float(s90["contact_sensor_angular_distance_deg"]) else 90.0
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
                "nearest_S0_S90_local_deg": near_0_90,
                # Positive means SENSOR_0's peak occurs earlier than SENSOR_90's.
                "S0_peak_lead_us": float(
                    (s90["peak_lag_s"] - s0["peak_lag_s"]) * 1e6
                ),
                "S0_crossing_lead_us": float(
                    (s90["first_5pct_peak_crossing_lag_s"]
                     - s0["first_5pct_peak_crossing_lag_s"])
                    * 1e6
                ),
            }
        )
    return pd.DataFrame(rows)


def plot_geometry(events: pd.DataFrame, summary: pd.DataFrame) -> None:
    """One ring diagram showing the four recurrent active-planet PR positions."""

    from matplotlib.patches import Circle

    fig, ax = plt.subplots(figsize=(6.5, 5.0))
    ax.set_aspect("equal")
    ax.set_xlim(-1.8, 1.8)
    ax.set_ylim(-1.38, 1.48)
    ax.axis("off")
    ax.add_patch(Circle((0, 0), 1.0, fill=False, lw=1.25, edgecolor="#313131"))
    ax.add_patch(Circle((0, 0), 0.88, fill=False, lw=0.65, edgecolor="#b7b7b7"))
    ax.plot([-1.07, 1.07], [0, 0], color="#dddddd", lw=0.55)
    ax.plot([0, 0], [-1.07, 1.07], color="#dddddd", lw=0.55)
    sensor_global: dict[float, float] = {}
    for sensor, group in events.groupby("sensor_angle_local_deg"):
        sensor_global[float(sensor)] = float(group.iloc[0]["sensor_angle_global_deg"])
    sensor_label_positions = {
        0.0: (0.0, 1.28, "center"),
        90.0: (1.23, 0.03, "left"),
        120.0: (1.16, -0.67, "left"),
        240.0: (-1.16, -0.67, "right"),
    }
    for sensor in SENSOR_ORDER:
        theta = np.deg2rad(sensor_global[sensor])
        x, y = 1.08 * np.cos(theta), 1.08 * np.sin(theta)
        ax.plot(x, y, marker="s", ms=5.6, color=SENSOR_COLORS[sensor], mec="white", mew=0.6)
        tx, ty, ha = sensor_label_positions[sensor]
        ax.text(
            tx,
            ty,
            f"S{int(sensor)} (FE {int(round(sensor_global[sensor])) % 360}°)",
            color=SENSOR_COLORS[sensor],
            fontsize=8.0,
            ha=ha,
            va="center",
        )

    for state, (letter, _target) in STATE_ORDER.items():
        representative = summary.loc[summary["state"].eq(state)].iloc[0]
        contact = float(representative["contact_angle_global_deg"])
        theta_contact = np.deg2rad(contact)
        x, y = np.cos(theta_contact), np.sin(theta_contact)
        ax.plot(
            [0, 0.96 * x],
            [0, 0.96 * y],
            color="#d62728",
            lw=0.55,
            ls="--",
            alpha=0.65,
            zorder=2,
        )
        ax.plot(
            0.96 * x,
            0.96 * y,
            marker="*",
            ms=10.5,
            color="#d62728",
            mec="white",
            mew=0.5,
            zorder=5,
        )
        ax.text(
            0.68 * x,
            0.68 * y,
            f"{letter}\n{contact:.1f}°",
            ha="center",
            va="center",
            fontsize=8.0,
            color="#b51f1f",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 0.8},
        )
    ax.annotate("FE angle +", xy=(0.26, 0.23), xytext=(0.36, 0.26), fontsize=7.2, color="#333333")
    fig.suptitle("Four recurring active-planet PR contact positions", fontsize=10, y=0.96)
    fig.text(0.5, 0.05, "Square: exported sensor patch; red star A–D: active-planet PR position at a sun-crack SP event.", ha="center", fontsize=7.2)
    fig.text(0.5, 0.025, "FE global angle increases counter-clockwise from the right; local sensor labels increase clockwise from S0.", ha="center", fontsize=7.2)
    fig.subplots_adjust(top=0.91, bottom=0.08)
    save_figure(fig, OUT / "Fig_FE_path_event_geometry", multi_panel=False)


def plot_sequence_table(summary: pd.DataFrame) -> None:
    """Export a compact event/planet/angle table without inventing a tooth ID."""

    # The source MAT contains a fault-event/planet index but no absolute tooth
    # number.  Make that limitation explicit instead of manufacturing a tooth
    # sequence from the event count.
    first = summary.head(16).copy()
    first["contact_angle_global_deg"] = first["contact_angle_global_deg"].map(lambda x: f"{x:.2f}")
    first["S0_peak_lead_us"] = first["S0_peak_lead_us"].map(lambda x: f"{x:+.0f}")
    first["nearest_sensor_local_deg"] = first["nearest_sensor_local_deg"].map(lambda x: f"S{int(x)}")
    first["strongest_sensor_local_deg"] = first["strongest_sensor_local_deg"].map(lambda x: f"S{int(x)}")
    table_rows = []
    for _, row in first.iterrows():
        table_rows.append(
            [
                f"{int(row['event_number'])}",
                f"P{int(row['active_planet'])}",
                f"{row['contact_angle_global_deg']}°",
                row["nearest_sensor_local_deg"],
                row["strongest_sensor_local_deg"],
                f"{row['S0_peak_lead_us']} μs",
            ]
        )

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.axis("off")
    columns = [
        "event",
        "active SP planet",
        "FE contact angle",
        "nearest",
        "max |Δy|",
        "S0 peak lead",
    ]
    table = ax.table(
        cellText=table_rows,
        colLabels=columns,
        cellLoc="center",
        colLoc="center",
        loc="upper center",
        bbox=[0.02, 0.11, 0.96, 0.82],
        colWidths=[0.08, 0.16, 0.19, 0.13, 0.13, 0.20],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7.0)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#b8b8b8")
        cell.set_linewidth(0.45)
        if row == 0:
            cell.set_facecolor("#e8eef4")
            cell.set_text_props(weight="bold")
        elif row % 4 == 0:
            cell.set_facecolor("#f6f6f6")
    fig.suptitle("Saved-source event sequence (first 16 events)", fontsize=10, y=0.965)
    fig.text(
        0.5,
        0.02,
        "Source records identify the active planet, not an absolute tooth number. No tooth ID is inferred.\n"
        "Positive S0 peak lead means SENSOR_0 reaches its modal-ringdown peak earlier.",
        ha="center",
        va="bottom",
        fontsize=7.0,
    )
    fig.subplots_adjust(top=0.91, bottom=0.08)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "Fig_FE_path_event_sequence.svg", bbox_inches="tight")
    fig.savefig(OUT / "Fig_FE_path_event_sequence.pdf", bbox_inches="tight")
    fig.savefig(OUT / "Fig_FE_path_event_sequence.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def plot_phase_lead(events: pd.DataFrame, signals: pd.DataFrame, summary: pd.DataFrame) -> None:
    """Two unshifted examples: near point can be strong but later at its peak."""

    time = signals["time_s"].to_numpy(float)
    y0 = signals["q0900um_minus_healthy_sensor_000deg_acceleration_ms2"].to_numpy(float)
    y90 = signals["q0900um_minus_healthy_sensor_090deg_acceleration_ms2"].to_numpy(float)
    selected = [3, 4, 8]
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 7.0), sharex=True, sharey=True)
    windows: list[tuple[np.ndarray, np.ndarray, np.ndarray, pd.Series, pd.Series]] = []
    ymax = 0.0
    for ev in selected:
        rows = events.loc[events["event_number"].eq(ev)]
        start = float(rows.iloc[0]["source_event_start_time_s"])
        keep = (time >= start) & (time <= start + 0.008)
        t_ms = (time[keep] - start) * 1000.0
        yy0, yy90 = y0[keep] * 1000.0, y90[keep] * 1000.0
        ymax = max(ymax, float(np.max(np.abs(np.r_[yy0, yy90]))))
        r0 = rows.loc[rows["sensor_angle_local_deg"].eq(0.0)].iloc[0]
        r90 = rows.loc[rows["sensor_angle_local_deg"].eq(90.0)].iloc[0]
        windows.append((t_ms, yy0, yy90, r0, r90))

    for ax, letter, ev, (t_ms, yy0, yy90, r0, r90) in zip(axes, ["a", "b", "c"], selected, windows):
        ax.plot(t_ms, yy0, color=SENSOR_COLORS[0.0], lw=0.85, label="SENSOR_0")
        ax.plot(t_ms, yy90, color=SENSOR_COLORS[90.0], lw=0.85, label="SENSOR_90")
        p0 = float(r0["peak_lag_s"]) * 1000.0
        p90 = float(r90["peak_lag_s"]) * 1000.0
        ax.axvline(p0, color=SENSOR_COLORS[0.0], ls=":", lw=0.9)
        ax.axvline(p90, color=SENSOR_COLORS[90.0], ls=":", lw=0.9)
        ax.axhline(0.0, color="#777777", lw=0.45)
        ax.grid(alpha=0.17, lw=0.4)
        angle = float(r0["contact_angle_global_deg"])
        near_0_90 = "S90" if angle > 300 else "S0"
        lead_us = (p0 - p90) * 1000.0
        if lead_us >= 0:
            lead_text = f"S90 peak {lead_us:.0f} μs before S0"
        else:
            lead_text = f"S0 peak {-lead_us:.0f} μs before S90"
        ax.set_title(
            f"({letter}) Event {ev}: FE contact {angle:.1f}°, nearer {near_0_90}; "
            + lead_text,
            fontsize=8.6,
            loc="left",
            pad=7,
        )
        ax.set_ylabel("Fault increment / 10⁻³ m s⁻²")
        ax.set_xlim(0.0, 8.0)
        ax.set_ylim(-1.30 * ymax, 1.30 * ymax)

    axes[-1].set_xlabel("Time from source-event start / ms")
    fig.suptitle("Unshifted 0°/90° sensor responses and modal-peak order", fontsize=10, y=0.98)
    fig.legend(
        handles=axes[0].lines[:2],
        labels=["SENSOR_0", "SENSOR_90"],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=2,
        frameon=False,
    )
    fig.text(
        0.5,
        0.02,
        "Dotted vertical lines mark each channel's event-window maximum |fault increment|; they do not mark the elastic-wave first arrival.",
        ha="center",
        va="bottom",
        fontsize=7.1,
    )
    fig.subplots_adjust(top=0.86, bottom=0.10, left=0.13, right=0.985, hspace=0.39)
    save_figure(fig, OUT / "Fig_FE_path_phase_lead_S0_S90", multi_panel=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    events, signals = load_data()
    summary = event_summary(events)
    summary.to_csv(OUT / "event_geometry_phase_summary.csv", index=False)
    plot_geometry(events, summary)
    plot_sequence_table(summary)
    plot_phase_lead(events, signals, summary)
    print(f"Wrote FE mapping evidence to {OUT}")
    print(
        "Counts: nearest=strongest %d/%d; 0/90 peak-lag near point earlier %d/%d; "
        "0/90 5%% crossing near point earlier %d/%d"
        % (
            int(summary["nearest_equals_strongest"].sum()),
            len(summary),
            int(sum(
                (
                    ((summary["nearest_S0_S90_local_deg"] == 0) & (summary["S0_peak_lead_us"] > 0))
                    | ((summary["nearest_S0_S90_local_deg"] == 90) & (summary["S0_peak_lead_us"] < 0))
                )
            )),
            len(summary),
            int(sum(
                (
                    ((summary["nearest_S0_S90_local_deg"] == 0) & (summary["S0_crossing_lead_us"] > 0))
                    | ((summary["nearest_S0_S90_local_deg"] == 90) & (summary["S0_crossing_lead_us"] < 0))
                )
            )),
            len(summary),
        )
    )


if __name__ == "__main__":
    main()
