"""Publication figures for the legacy Chapter-3 18+7 DOF route.

The script reads the *saved* paper-v3 three-channel response and performs only
post-processing.  It does not call the FE ring transfer path and it does not
alter the source dynamics.  Outputs are intentionally written to a separate
results directory so that historical figures remain untouched.

Example
-------
python scripts/plot_ch3_old18plus7_figures.py \
  --signals "D:/.../results/paper_v3_three_channel_signals_nominal.csv" \
  --out "D:/.../results/chapter3_old18plus7_figures"
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


try:
    from audit_panel_alignment import require_matplotlib_panel_alignment
except ImportError:  # pragma: no cover - the script remains usable without QA helper
    require_matplotlib_panel_alignment = None


FS_DEFAULT = 51200.0
FM_HZ = 168.0
SENSOR_LABELS = ("0°", "120°", "240°")
SENSOR_KEYS = ("000", "120", "240")
CASES = ("healthy", "sun25", "sun50")
CASE_LABELS = ("Healthy", "Sun crack 25%", "Sun crack 50%")
CASE_COLORS = ("#555555", "#1769AA", "#C43C39")
SENSOR_COLORS = ("#1B5EAA", "#D97706", "#238B45")


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            # Keep a publication-safe sans fallback for systems without Times.
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "mathtext.fontset": "stix",
            "font.size": 8.0,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9.0,
            "xtick.labelsize": 7.2,
            "ytick.labelsize": 7.2,
            "legend.fontsize": 7.2,
            "axes.linewidth": 0.75,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "xtick.major.size": 3.2,
            "ytick.major.size": 3.2,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.11,
        1.05,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontweight="bold",
        fontsize=9,
    )


def save_bundle(fig: plt.Figure, stem: Path, *, multipanel: bool = False) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    qa_dir = stem.parent / "qa"
    qa_dir.mkdir(parents=True, exist_ok=True)
    if multipanel and require_matplotlib_panel_alignment is not None:
        require_matplotlib_panel_alignment(
            fig,
            json_out=str(qa_dir / f"{stem.name}.alignment.json"),
            overlay_svg=str(qa_dir / f"{stem.name}.alignment.svg"),
            tolerance_pt=1.5,
            gutter_tolerance_pt=1.5,
            require_panel_labels=True,
            strict=True,
        )
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def read_signals(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    frame = pd.read_csv(path)
    if "time_s" not in frame:
        raise ValueError("signal CSV must contain time_s")
    time = frame["time_s"].to_numpy(float)
    signals: dict[str, np.ndarray] = {}
    for case in CASES:
        for sensor in SENSOR_KEYS:
            name = f"{case}_sensor_{sensor}deg_ms2"
            if name not in frame:
                raise ValueError(f"missing column: {name}")
            signals[name] = frame[name].to_numpy(float)
    if len(time) < 16 or not np.all(np.isfinite(time)):
        raise ValueError("invalid time vector")
    dt = np.diff(time)
    if not np.allclose(dt, np.median(dt), rtol=2e-6, atol=1e-12):
        raise ValueError("time vector is not uniformly sampled")
    return time, signals


def one_sided_fft(x: np.ndarray, fs: float) -> tuple[np.ndarray, np.ndarray]:
    y = np.asarray(x, float) - np.mean(x)
    win = np.hanning(len(y))
    spec = 2.0 * np.abs(np.fft.rfft(y * win)) / np.sum(win)
    freq = np.fft.rfftfreq(len(y), 1.0 / fs)
    if spec.size:
        spec[0] *= 0.5
    return freq, spec


def compute_spectra(time: np.ndarray, signals: dict[str, np.ndarray]) -> pd.DataFrame:
    fs = 1.0 / np.median(np.diff(time))
    rows: list[pd.DataFrame] = []
    for case, label in zip(CASES, CASE_LABELS):
        for sensor in SENSOR_KEYS:
            key = f"{case}_sensor_{sensor}deg_ms2"
            freq, amp = one_sided_fft(signals[key], fs)
            rows.append(
                pd.DataFrame(
                    {
                        "case": case,
                        "case_label": label,
                        "sensor_deg": int(sensor),
                        "frequency_Hz": freq,
                        "amplitude": amp,
                    }
                )
            )
    return pd.concat(rows, ignore_index=True)


def plot_time(time: np.ndarray, signals: dict[str, np.ndarray], out: Path) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 6.3), sharex=True,
                             gridspec_kw={"hspace": 0.22})
    for i, (case, label) in enumerate(zip(CASES, CASE_LABELS)):
        ax = axes[i]
        for sensor, colour, sensor_label in zip(SENSOR_KEYS, SENSOR_COLORS, SENSOR_LABELS):
            key = f"{case}_sensor_{sensor}deg_ms2"
            ax.plot(time, signals[key], color=colour, lw=0.55, label=sensor_label)
        ax.set_ylabel("Amplitude")
        ax.set_title(label, loc="left", pad=3)
        ax.axhline(0.0, color="#B8B8B8", lw=0.45, zorder=0)
        ax.grid(axis="y", color="#E6E6E6", lw=0.45)
        ax.set_xlim(float(time[0]), float(time[-1]))
        panel_label(ax, "abc"[i])
        if i == 0:
            ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02),
                      ncol=3, handlelength=2.0, borderaxespad=0.0)
    axes[-1].set_xlabel("Time / s")
    fig.suptitle("System vibration response", y=0.995, fontsize=10)
    fig.subplots_adjust(left=0.105, right=0.985, top=0.94, bottom=0.085)
    save_bundle(fig, out / "Fig3_old18plus7_system_response_time", multipanel=True)


def plot_spectrum(spec: pd.DataFrame, out: Path, xmax: float = 600.0) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 6.3), sharex=True,
                             gridspec_kw={"hspace": 0.22})
    for i, (case, label) in enumerate(zip(CASES, CASE_LABELS)):
        ax = axes[i]
        for sensor, colour, sensor_label in zip(SENSOR_KEYS, SENSOR_COLORS, SENSOR_LABELS):
            m = (spec.case == case) & (spec.sensor_deg == int(sensor)) & (spec.frequency_Hz <= xmax)
            ax.plot(spec.loc[m, "frequency_Hz"], spec.loc[m, "amplitude"],
                    color=colour, lw=0.75, label=sensor_label)
        for k, colour, text in ((1, "#6B7280", "fₘ"), (2, "#9CA3AF", "2fₘ"), (3, "#D1D5DB", "3fₘ")):
            f0 = k * FM_HZ
            if f0 <= xmax:
                ax.axvline(f0, color=colour, lw=0.6, ls=(0, (2, 2)), zorder=0)
        ax.set_ylabel("Amplitude")
        ax.set_title(label, loc="left", pad=3)
        ax.set_xlim(0, xmax)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", color="#E6E6E6", lw=0.45)
        panel_label(ax, "abc"[i])
        if i == 0:
            ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02),
                      ncol=3, handlelength=2.0, borderaxespad=0.0)
    axes[-1].set_xlabel("Frequency / Hz")
    fig.suptitle("System vibration response spectra", y=0.995, fontsize=10)
    fig.subplots_adjust(left=0.105, right=0.985, top=0.94, bottom=0.085)
    save_bundle(fig, out / "Fig3_old18plus7_system_response_spectrum_0_600Hz", multipanel=True)


def _downsample(x: np.ndarray, y: np.ndarray, max_points: int = 6000) -> tuple[np.ndarray, np.ndarray]:
    if len(x) <= max_points:
        return x, y
    idx = np.linspace(0, len(x) - 1, max_points, dtype=int)
    return x[idx], y[idx]


def plot_waterfalls(time: np.ndarray, signals: dict[str, np.ndarray], spec: pd.DataFrame, out: Path) -> None:
    severity = np.array([0.0, 25.0, 50.0])
    for sensor, sensor_label in zip(SENSOR_KEYS, SENSOR_LABELS):
        # Time waterfall
        fig = plt.figure(figsize=(7.2, 4.8), facecolor="white")
        ax = fig.add_subplot(111, projection="3d")
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.fill = False
            axis.pane.set_edgecolor("#D0D0D0")
        for case, sev, colour in zip(CASES, severity, CASE_COLORS):
            x, y = _downsample(time, signals[f"{case}_sensor_{sensor}deg_ms2"])
            ax.plot(x, np.full_like(x, sev), y, color=colour, lw=0.55)
        ax.set_xlabel("Time / s", labelpad=6)
        ax.set_ylabel("Sun crack / %", labelpad=12)
        # Matplotlib's 3-D z-label can be clipped by bbox_inches='tight'.
        # Put the editable label in figure coordinates so it remains visible
        # in PNG, PDF and SVG exports.
        ax.set_zlabel("")
        ax.set_yticks(severity)
        ax.set_yticklabels(["0", "25", "50"])
        ax.view_init(elev=18, azim=-68)
        ax.set_title(f"Time response — sensor {sensor_label}", pad=10)
        fig.subplots_adjust(left=0.03, right=0.86, bottom=0.07, top=0.88)
        fig.text(0.88, 0.50, "Amplitude", rotation=90, rotation_mode="anchor",
                 ha="center", va="center", fontsize=8.5)
        save_bundle(fig, out / f"Fig3_old18plus7_time_waterfall_sensor_{sensor}deg")

        # Spectrum waterfall (0–600 Hz; includes 1fm, 2fm and 3fm)
        fig = plt.figure(figsize=(7.2, 4.8), facecolor="white")
        ax = fig.add_subplot(111, projection="3d")
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.fill = False
            axis.pane.set_edgecolor("#D0D0D0")
        for case, sev, colour in zip(CASES, severity, CASE_COLORS):
            m = (spec.case == case) & (spec.sensor_deg == int(sensor)) & (spec.frequency_Hz <= 600.0)
            ax.plot(spec.loc[m, "frequency_Hz"], np.full(m.sum(), sev),
                    spec.loc[m, "amplitude"], color=colour, lw=0.65)
        for k in (1, 2, 3):
            ax.plot([k * FM_HZ, k * FM_HZ], [0, 50], [0, 0],
                    color="#B7B7B7", lw=0.45, ls=(0, (2, 2)))
        ax.set_xlabel("Frequency / Hz", labelpad=6)
        ax.set_ylabel("Sun crack / %", labelpad=12)
        ax.set_zlabel("")
        ax.set_yticks(severity)
        ax.set_yticklabels(["0", "25", "50"])
        ax.set_xlim(0, 600)
        ax.view_init(elev=18, azim=-68)
        ax.set_title(f"Response spectrum — sensor {sensor_label}", pad=10)
        fig.subplots_adjust(left=0.03, right=0.86, bottom=0.07, top=0.88)
        fig.text(0.88, 0.50, "Amplitude", rotation=90, rotation_mode="anchor",
                 ha="center", va="center", fontsize=8.5)
        save_bundle(fig, out / f"Fig3_old18plus7_spectrum_waterfall_sensor_{sensor}deg")


def write_metrics(time: np.ndarray, signals: dict[str, np.ndarray], spec: pd.DataFrame, out: Path) -> None:
    rows = []
    for case, label in zip(CASES, CASE_LABELS):
        for sensor in SENSOR_KEYS:
            x = signals[f"{case}_sensor_{sensor}deg_ms2"]
            rows.append({
                "case": case,
                "case_label": label,
                "sensor_deg": int(sensor),
                "mean": float(np.mean(x)),
                "rms": float(np.sqrt(np.mean(x**2))),
                "std": float(np.std(x)),
                "min": float(np.min(x)),
                "max": float(np.max(x)),
                "peak_to_peak": float(np.ptp(x)),
            })
    pd.DataFrame(rows).to_csv(out / "system_response_metrics.csv", index=False)
    spec.to_csv(out / "system_response_spectra_full.csv", index=False)
    spec[spec.frequency_Hz <= 600.0].to_csv(out / "system_response_spectra_0_600Hz.csv", index=False)


def write_readme(source: Path, out: Path, time: np.ndarray) -> None:
    fs = 1.0 / np.median(np.diff(time))
    text = f"""Chapter 3 old-route figure bundle

Source CSV: {source}
Model route: legacy paper-v3 18 rigid DOF + 7 effective analytical ring modes (25 total coordinates).
Sensors: fixed radial channels at 0, 120 and 240 degrees.
This bundle does not use the independent ANSYS FE one-way transfer-path route.

Sampling rate: {fs:.12g} Hz
Record interval: {time[0]:.9g} to {time[-1]:.9g} s
Mesh frequency used for reference markers: {FM_HZ:g} Hz

The time traces are the saved model acceleration responses.  The plotted y-axis
is intentionally displayed as `Amplitude`; source values remain unchanged.
Spectra use mean removal, a Hann window, and one-sided amplitude scaling.  No
extra fault pulses, sidebands, smoothing, or frequency lines are inserted.

Render-time alignment and collision QA artifacts are stored in the `qa/`
subdirectory beside the figure files.

Operating inputs in the source model are provisional engineering settings (20 N m
input torque, 100 N m ideal carrier load and 5 um TE); they are not measured data.
"""
    (out / "README.txt").write_text(text, encoding="utf-8")


def plot_path_truth(path_csv: Path, out: Path) -> None:
    """Plot the saved legacy path-truth FRF, without interpreting it as FE data."""
    frame = pd.read_csv(path_csv)
    frequencies = sorted(frame["frequency_Hz"].unique())
    if not frequencies:
        return
    # The source grid used for this diagnostic is 1.7--1.848 kHz; retain the
    # actual grid values rather than relabelling them as the 168-Hz mesh line.
    chosen = [frequencies[0], frequencies[-1]] if len(frequencies) > 1 else [frequencies[0]]
    fig, axes = plt.subplots(2, 2 if len(chosen) > 1 else 1,
                             figsize=(7.2, 4.9), squeeze=False,
                             gridspec_kw={"wspace": 0.28, "hspace": 0.36})
    for j, freq in enumerate(chosen):
        sub = frame[np.isclose(frame.frequency_Hz, freq)]
        ax_amp = axes[0, j]
        ax_phase = axes[1, j]
        for sensor, colour, label in zip(
            ("sensor_000deg", "sensor_120deg", "sensor_240deg"),
            SENSOR_COLORS,
            SENSOR_LABELS,
        ):
            s = sub[sub.sensor_label == sensor].sort_values("carrier_angle_deg")
            ax_amp.plot(s.carrier_angle_deg, s.amplitude, color=colour, lw=0.9,
                        marker="o", ms=2.0, label=label)
            ax_phase.plot(s.carrier_angle_deg, s.phase_difference_deg, color=colour,
                          lw=0.9, marker="o", ms=2.0, label=label)
        ax_amp.set_title(f"{freq:.0f} Hz", loc="left", pad=3)
        ax_amp.set_ylabel("Path amplitude")
        ax_amp.set_xlim(0, 360)
        ax_amp.set_xticks([0, 90, 180, 270, 360])
        ax_amp.grid(axis="y", color="#E6E6E6", lw=0.45)
        ax_phase.set_xlabel("Carrier angle / deg")
        ax_phase.set_ylabel("Phase difference / deg")
        ax_phase.set_xlim(0, 360)
        ax_phase.set_xticks([0, 90, 180, 270, 360])
        ax_phase.axhline(0, color="#B8B8B8", lw=0.45)
        ax_phase.grid(axis="y", color="#E6E6E6", lw=0.45)
        panel_label(ax_amp, "abcd"[2*j])
        panel_label(ax_phase, "abcd"[2*j+1])
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.955),
               ncol=3, handlelength=2.0, frameon=False)
    fig.suptitle("Transfer-path amplitude and phase response", y=0.995, fontsize=10)
    fig.subplots_adjust(left=0.09, right=0.985, top=0.86, bottom=0.12)
    save_bundle(fig, out / "Fig3_old18plus7_path_amplitude_phase", multipanel=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signals", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    configure_style()
    time, signals = read_signals(args.signals)
    spec = compute_spectra(time, signals)
    args.out.mkdir(parents=True, exist_ok=True)
    plot_time(time, signals, args.out)
    plot_spectrum(spec, args.out)
    plot_waterfalls(time, signals, spec, args.out)
    path_csv = args.out / "path_truth_sensor_event.csv"
    if path_csv.exists():
        plot_path_truth(path_csv, args.out)
    write_metrics(time, signals, spec, args.out)
    write_readme(args.signals, args.out, time)
    meta = {
        "source": str(args.signals),
        "route": "legacy 18 rigid DOF + 7 effective ring modes",
        "sensors_deg": [0, 120, 240],
        "mesh_frequency_Hz": FM_HZ,
        "fft": "mean removed + Hann + one-sided amplitude",
    }
    (args.out / "figure_manifest.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"Saved Chapter-3 figures to {args.out}")


if __name__ == "__main__":
    main()
