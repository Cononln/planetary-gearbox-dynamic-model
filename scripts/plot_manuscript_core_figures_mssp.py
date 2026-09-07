"""Create the workflow and three core MSSP result figures from frozen source data."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "paper_figures_mssp" / "source_data"
OUT = ROOT / "results" / "paper_figures_mssp"

NAVY = "#264653"
BLUE = "#3A6EA5"
TEAL = "#2A9D8F"
GOLD = "#E9C46A"
RED = "#D1495B"
LIGHT_BLUE = "#A8C4E0"
MID_GREY = "#747474"
LIGHT_GREY = "#E7E7E7"
INK = "#222222"

METHOD_COLORS = {
    "single_channel_H1": "#A6A6A6",
    "multiharmonic_no_path": LIGHT_BLUE,
    "unconstrained_joint": GOLD,
    "hard_three_planet_joint": TEAL,
    "adaptive_three_planet_joint": RED,
}

METHOD_LABELS = {
    "single_channel_H1": "Single H1",
    "multiharmonic_no_path": "Multi-H",
    "unconstrained_joint": "Free path",
    "hard_three_planet_joint": "Strict 3k",
    "adaptive_three_planet_joint": "Adaptive 3k",
}

DIAG_COLORS = [INK, BLUE, RED]
DIAG_LABELS = ["Single sensor", "Naive fusion", "Phase-aligned fusion"]


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 7.2,
            "axes.labelsize": 7.5,
            "axes.titlesize": 8.0,
            "xtick.labelsize": 6.7,
            "ytick.labelsize": 6.7,
            "legend.fontsize": 6.6,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.06) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontweight="bold", fontsize=8.5,
            va="top", ha="left")


def save_bundle(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def rounded_box(ax: plt.Axes, xy: tuple[float, float], wh: tuple[float, float],
                title: str, lines: list[str], edge: str, fill: str) -> None:
    x, y = xy
    w, h = wh
    box = FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.025",
        linewidth=1.0, edgecolor=edge, facecolor=fill,
    )
    ax.add_patch(box)
    ax.text(x + 0.04 * w, y + 0.77 * h, title, color=edge, fontweight="bold",
            fontsize=7.4, ha="left", va="center")
    ax.text(x + 0.04 * w, y + 0.47 * h, "\n".join(lines), color=INK,
            fontsize=6.5, ha="left", va="center", linespacing=1.22)


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float],
          color: str = MID_GREY, rad: float = 0.0) -> None:
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=8,
                                 linewidth=0.9, color=color,
                                 connectionstyle=f"arc3,rad={rad}"))


def plot_workflow() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.15))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Physical observation panel
    ax.text(0.02, 0.96, "a", fontweight="bold", fontsize=8.5, va="top")
    ax.text(0.055, 0.96, "Moving-path observation", fontweight="bold", va="top")
    cx, cy = 0.12, 0.58
    ax.add_patch(Circle((cx, cy), 0.105, fill=False, lw=2.0, ec=MID_GREY))
    ax.add_patch(Circle((cx, cy), 0.025, fc=GOLD, ec=INK, lw=0.7))
    for ang in np.deg2rad([20, 140, 260]):
        px, py = cx + 0.064 * np.cos(ang), cy + 0.064 * np.sin(ang)
        ax.add_patch(Circle((px, py), 0.018, fc=LIGHT_BLUE, ec=BLUE, lw=0.8))
    sensors = [(cx + 0.125, cy, "S1"), (cx, cy + 0.125, "S2"),
               (cx - 0.108, cy - 0.07, "S3")]
    for sx, sy, name in sensors:
        ax.plot(sx, sy, marker="^", ms=6, color=RED, mec=INK, mew=0.5)
        ax.text(sx + 0.012, sy + 0.008, name, fontsize=6.2)
    px, py = cx + 0.064 * np.cos(np.deg2rad(20)), cy + 0.064 * np.sin(np.deg2rad(20))
    arrow(ax, (px, py), (sensors[0][0] - 0.01, sensors[0][1]), BLUE, 0.12)
    arrow(ax, (px, py), (sensors[1][0], sensors[1][1] - 0.01), RED, -0.12)
    ax.text(0.025, 0.31, "Common speed phase +\nsensor-dependent path phase",
            fontsize=6.5, ha="left", va="top", color=INK)
    ax.text(0.025, 0.14, "Input: synchronized vibration only", fontsize=6.3,
            color=NAVY, fontweight="bold")
    ax.text(0.025, 0.075, "Encoder hidden during estimation", fontsize=6.1,
            color=MID_GREY)

    # Method modules
    rounded_box(ax, (0.27, 0.58), (0.20, 0.28), "b  Common phase",
                ["Multi-sensor, multi-harmonic", "phase increments", r"$\hat{\phi}_m \rightarrow \hat{\theta}_c=\hat{\phi}_m/Z_r$"],
                NAVY, "#EEF3F5")
    rounded_box(ax, (0.52, 0.58), (0.20, 0.28), "c  Differential path",
                ["Pairwise cross-phasors", "Three-planet spatial orders", "Zero-mean sensor graph"],
                TEAL, "#EEF7F5")
    rounded_box(ax, (0.77, 0.58), (0.20, 0.28), "d  Topology safeguard",
                ["Odd/even carrier-cycle holdout", "Complete graph / reliable pair", "No-path return if unsupported"],
                RED, "#FAF0F1")
    arrow(ax, (0.225, 0.70), (0.27, 0.70))
    arrow(ax, (0.47, 0.72), (0.52, 0.72))
    arrow(ax, (0.72, 0.72), (0.77, 0.72))

    # Output row
    rounded_box(ax, (0.29, 0.12), (0.25, 0.25), "e  Unit-modulus correction",
                [r"$\tilde z_{s,h}=z_{s,h}\exp(-j\hat\psi_{s,h})$", r"$|\tilde z_{s,h}|=|z_{s,h}|$", "Relative phase only; amplitude retained"],
                BLUE, "#F0F4FA")
    rounded_box(ax, (0.61, 0.12), (0.30, 0.25), "f  Phase-aware outputs",
                ["Coherent common mode", "Differential residuals", "Envelope spectra and TPSVD"],
                NAVY, "#F3F3F3")
    arrow(ax, (0.82, 0.58), (0.52, 0.37), RED, 0.0)
    arrow(ax, (0.54, 0.245), (0.61, 0.245))

    save_bundle(fig, "method_workflow_mssp")


def break_wrapped_curve(angle_deg: np.ndarray, phase_deg: np.ndarray,
                        jump: float = 110.0) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(angle_deg, float).copy()
    y = np.asarray(phase_deg, float).copy()
    bad = np.abs(np.diff(y)) > jump
    x2, y2 = [x[0]], [y[0]]
    for i in range(1, len(x)):
        if bad[i - 1]:
            x2.append(np.nan)
            y2.append(np.nan)
        x2.append(x[i])
        y2.append(y[i])
    return np.asarray(x2), np.asarray(y2)


def plot_dynamic_path() -> None:
    data = pd.read_csv(SOURCE / "dynamic_path_field.csv")
    metric = pd.read_csv(ROOT / "results" / "loaded_path_am_pm_metrics.csv")
    valid = pd.read_csv(ROOT / "results" / "loaded_path_am_pm_validation.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.65), sharex=True,
                             gridspec_kw={"wspace": 0.25, "hspace": 0.34})
    freqs = [168, 1848]
    for row, freq in enumerate(freqs):
        ax_a, ax_p = axes[row]
        a0 = data[f"amp_{freq}_sensor0"].to_numpy()
        a90 = data[f"amp_{freq}_sensor90"].to_numpy()
        scale = max(np.nanmax(a0), np.nanmax(a90))
        angle = data["carrier_angle_deg"].to_numpy()
        ax_a.plot(angle, a0 / scale, color=NAVY, lw=1.3, label="Sensor at 0°")
        ax_a.plot(angle, a90 / scale, color=RED, lw=1.3, label="Sensor at 90°")
        ax_a.set_ylabel("Normalized path magnitude")
        ax_a.set_ylim(-0.02, 1.07)
        ax_a.set_yticks([0, 0.5, 1.0])
        ax_a.set_title(f"{freq} Hz: path magnitude", loc="left", pad=4)
        if row == 0:
            ax_a.legend(loc="upper right", ncol=1, handlelength=2.2)

        xp, yp = break_wrapped_curve(angle, data[f"phase_{freq}_deg"].to_numpy())
        ax_p.plot(xp, yp, color=INK, lw=1.3)
        ax_p.axhline(0, color=LIGHT_GREY, lw=0.7, zorder=0)
        ax_p.set_ylim(-190, 190)
        ax_p.set_yticks([-180, -90, 0, 90, 180])
        ax_p.set_ylabel("S2–S1 phase (deg)")
        ax_p.set_title(f"{freq} Hz: differential path phase", loc="left", pad=4)
        span = metric.loc[metric.frequency_hz.eq(freq), "dual_phase_span_deg"].iloc[0]
        corr = valid.loc[valid.frequency_hz.eq(freq),
                         "path_simulation_sideband_correlation"].iloc[0]
        ax_p.text(0.03, 0.08, f"phase span = {span:.1f}°\nsideband r = {corr:.4f}",
                  transform=ax_p.transAxes, ha="left", va="bottom", color=NAVY,
                  fontsize=6.5)

        for ax in (ax_a, ax_p):
            ax.set_xlim(0, 360)
            ax.set_xticks([0, 90, 180, 270, 360])
            ax.tick_params(direction="out")
            if row == 1:
                ax.set_xlabel("Carrier angle (deg)")

    for label, ax in zip("abcd", axes.flat):
        panel_label(ax, label)
    save_bundle(fig, "dynamic_path_am_pm_mssp")


def plot_encoder_hidden_validation() -> None:
    freq = pd.read_csv(SOURCE / "encoder_hidden_frequency.csv")
    path = pd.read_csv(SOURCE / "encoder_hidden_path.csv")
    metrics = pd.read_csv(SOURCE / "encoder_hidden_metrics.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.75),
                             gridspec_kw={"wspace": 0.28, "hspace": 0.38})

    ax = axes[0, 0]
    ax.plot(freq.time_s, freq.encoder_mesh_frequency_hz, color=INK, lw=1.0,
            label="Hidden encoder")
    ax.plot(freq.time_s, freq.blind_mesh_frequency_hz, color=RED, lw=1.0,
            label="Vibration-only estimate")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Mesh frequency (Hz)")
    ax.set_title("Blind common-phase tracking", loc="left", pad=4)
    ax.legend(loc="upper left")
    adaptive = metrics.loc[metrics.method.eq("adaptive_three_planet_joint")].iloc[0]
    ax.text(0.98, 0.08, f"carrier-speed RMSE\n{adaptive.carrier_speed_rmse_rpm:.3f} r/min",
            transform=ax.transAxes, ha="right", va="bottom", color=NAVY)

    ax = axes[0, 1]
    xo, yo = break_wrapped_curve(path.carrier_angle_deg.to_numpy(),
                                 path.encoder_oracle_path_deg.to_numpy())
    xb, yb = break_wrapped_curve(path.carrier_angle_deg.to_numpy(),
                                 path.blind_path_deg.to_numpy())
    ax.plot(xo, yo, color=INK, lw=1.25, label="Encoder-angle oracle")
    ax.plot(xb, yb, color=RED, lw=1.0, label="Blind estimate")
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_ylim(-190, 190)
    ax.set_yticks([-180, -90, 0, 90, 180])
    ax.set_xlabel("Carrier angle (deg)")
    ax.set_ylabel("S2–S1 path phase (deg)")
    ax.set_title("Differential-path external score", loc="left", pad=4)
    ax.legend(loc="lower left")
    ax.text(0.98, 0.08, f"path RMSE = {adaptive.path_phase_vs_encoder_rmse_deg:.2f}°",
            transform=ax.transAxes, ha="right", va="bottom", color=NAVY)

    ax = axes[1, 0]
    for _, row in metrics.iterrows():
        ax.scatter(row.dual_phase_residual_deg, row.dual_phase_coherence, s=34,
                   color=METHOD_COLORS[row.method], edgecolor="white", linewidth=0.5,
                   zorder=3)
    no_path = metrics.loc[metrics.method.eq("multiharmonic_no_path")].iloc[0]
    free = metrics.loc[metrics.method.eq("unconstrained_joint")].iloc[0]
    strict = metrics.loc[metrics.method.eq("hard_three_planet_joint")].iloc[0]
    adaptive = metrics.loc[metrics.method.eq("adaptive_three_planet_joint")].iloc[0]
    ax.text(no_path.dual_phase_residual_deg - 8.0, no_path.dual_phase_coherence + 0.010,
            "H1 / Multi-H (no path)", fontsize=6.2, color=INK)
    ax.text(strict.dual_phase_residual_deg - 12.5, strict.dual_phase_coherence + 0.012,
            "Strict 3k", fontsize=6.2, color=INK)
    ax.text(free.dual_phase_residual_deg + 0.8, free.dual_phase_coherence + 0.013,
            "Free path", fontsize=6.2, color=INK)
    ax.text(adaptive.dual_phase_residual_deg + 0.8, adaptive.dual_phase_coherence - 0.020,
            "Adaptive 3k", fontsize=6.2, color=INK)
    ax.set_xlabel("Inter-channel phase residual (deg; lower is better)")
    ax.set_ylabel("Inter-channel coherence (higher is better)")
    ax.set_title("Alignment trade-off", loc="left", pad=4)
    ax.set_xlim(44, 74)
    ax.set_ylim(0.53, 0.83)

    ax = axes[1, 1]
    joint = metrics[metrics.method.isin(
        ["unconstrained_joint", "hard_three_planet_joint", "adaptive_three_planet_joint"]
    )].copy()
    x = np.arange(len(joint))
    bars = ax.bar(x, joint.path_phase_vs_encoder_rmse_deg,
                  color=[METHOD_COLORS[m] for m in joint.method], width=0.62,
                  edgecolor="white", linewidth=0.7)
    ax.set_xticks(x, [METHOD_LABELS[m] for m in joint.method])
    ax.set_ylabel("Path-phase RMSE (deg)")
    ax.set_title("Encoder-scored path recovery", loc="left", pad=4)
    ax.set_ylim(0, 55)
    for bar, value in zip(bars, joint.path_phase_vs_encoder_rmse_deg):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 1.2, f"{value:.1f}",
                ha="center", va="bottom", fontsize=6.5)
    ax.text(0.98, 0.93, "Encoder used for scoring only", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.2, color=MID_GREY)

    for label, ax in zip("abcd", axes.flat):
        panel_label(ax, label)
        ax.tick_params(direction="out")
    save_bundle(fig, "encoder_hidden_validation_mssp")


def read_cycle_matrix(index: int) -> np.ndarray:
    frame = pd.read_csv(SOURCE / f"sun_fault_cycle_matrix_{index}.csv")
    return frame.drop(columns=["cycle_index"]).to_numpy(float)


def cumulative_energy(values: np.ndarray) -> np.ndarray:
    values = np.nan_to_num(np.asarray(values, float), nan=0.0)
    power = values ** 2
    return np.cumsum(power) / np.sum(power)


def plot_sun_fault_diagnosis() -> None:
    spectrum = pd.read_csv(SOURCE / "sun_fault_envelope_spectrum.csv")
    metrics = pd.read_csv(SOURCE / "sun_fault_metrics.csv")
    sv = pd.read_csv(SOURCE / "sun_fault_singular_values.csv")
    patterns = pd.read_csv(SOURCE / "sun_fault_tpsvd_pattern.csv")
    naive = read_cycle_matrix(2)
    aligned = read_cycle_matrix(3)
    angle = pd.read_csv(SOURCE / "sun_fault_cycle_angle_deg.csv", header=None).iloc[:, 0].to_numpy()

    fig = plt.figure(figsize=(7.2, 5.15))
    gs = fig.add_gridspec(2, 3, wspace=0.34, hspace=0.40)
    axes = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(3)]

    ax = axes[0]
    cols = ["single_sensor", "naive_fusion", "phase_aligned_fusion"]
    for col, color, label in zip(cols, DIAG_COLORS, DIAG_LABELS):
        ax.plot(spectrum.frequency_hz, spectrum[col], color=color, lw=0.9, label=label)
    f_sf = 23.83
    for h in range(1, 5):
        ax.axvline(h * f_sf, color=MID_GREY, ls=(0, (2, 2)), lw=0.65, zorder=0)
    ax.set_xlim(0, 105)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Envelope frequency (Hz)")
    ax.set_ylabel("Amplitude")
    ax.set_title("Sun-fault family in the envelope", loc="left", pad=4)
    ax.legend(loc="upper right", handlelength=2.2)

    ax = axes[1]
    harmonics = np.arange(1, 5)
    width = 0.23
    for k, (row, color, label) in enumerate(zip(metrics.itertuples(), DIAG_COLORS, DIAG_LABELS)):
        values = np.array([row.fault_1x_amplitude, row.fault_2x_amplitude,
                           row.fault_3x_amplitude, row.fault_4x_amplitude])
        ax.bar(harmonics + (k - 1) * width, values, width=width, color=color,
               label=label, edgecolor="white", linewidth=0.4)
    ax.set_xticks(harmonics, [r"$1f_{sf}$", r"$2f_{sf}$", r"$3f_{sf}$", r"$4f_{sf}$"])
    ax.set_ylabel("Envelope amplitude")
    ax.set_title("Fault-line amplitudes", loc="left", pad=4)
    ax.text(0.98, 0.96, "one 50-s record; no error bars", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.1, color=MID_GREY)

    vmax = np.nanpercentile(np.abs(np.r_[naive.ravel(), aligned.ravel()]), 99.2)
    heat_axes = [axes[3], axes[4]]
    for ax, matrix, title in zip(heat_axes, [naive, aligned],
                                 ["Naive-fusion carrier cycles", "Aligned-fusion carrier cycles"]):
        im = ax.imshow(matrix, aspect="auto", origin="lower", cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax, extent=[angle.min(), angle.max(), 1, matrix.shape[0]])
        ax.set_xlabel("Carrier angle (deg)")
        ax.set_ylabel("Cycle index")
        ax.set_title(title, loc="left", pad=4)
        ax.set_xticks([0, 90, 180, 270, 360])
    cax = inset_axes(heat_axes[1], width="3.4%", height="42%", loc="upper right",
                     borderpad=0.8)
    cbar = fig.colorbar(im, cax=cax)
    cax.set_title("ampl.", fontsize=5.8, pad=2)
    cbar.ax.tick_params(labelsize=6.1)

    ax = axes[2]
    for col, color, label in zip(cols, DIAG_COLORS, DIAG_LABELS):
        vals = sv[col].dropna().to_numpy()
        rank = np.arange(1, min(12, len(vals)) + 1)
        ax.plot(rank, cumulative_energy(vals)[: len(rank)], color=color, lw=1.1,
                marker="o", ms=2.8, label=label)
    ax.set_xlim(1, 12)
    ax.set_ylim(0.68, 0.99)
    ax.set_xlabel("Retained singular modes")
    ax.set_ylabel("Cumulative TPSVD energy")
    ax.set_title("Periodic-mode compactness", loc="left", pad=4)
    ax.legend(loc="lower right", handlelength=2.0)

    ax = axes[5]
    n = len(patterns)
    order = np.fft.rfftfreq(n, d=1 / n)
    keep = order <= 40
    for col, color, label in zip(cols, DIAG_COLORS, DIAG_LABELS):
        y = patterns[col].to_numpy(float)
        y = y - np.nanmean(y)
        amp = 2.0 * np.abs(np.fft.rfft(y)) / n
        ax.plot(order[keep], amp[keep], color=color, lw=1.0, label=label)
    for o in [12, 24]:
        ax.axvline(o, color=MID_GREY, ls=(0, (2, 2)), lw=0.65, zorder=0)
    ax.text(12.5, 0.96, r"$f_{sf}/f_c$", transform=ax.get_xaxis_transform(),
            ha="left", va="top", fontsize=6.2)
    ax.set_xlim(0, 40)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Carrier order")
    ax.set_ylabel("Rank-1 pattern amplitude")
    ax.set_title("TPSVD fault-order pattern", loc="left", pad=4)

    for label, ax in zip("abcdef", axes):
        panel_label(ax, label)
        ax.tick_params(direction="out")
    save_bundle(fig, "sun_fault_diagnosis_tpsvd_mssp")


def main() -> None:
    configure_style()
    plot_workflow()
    plot_dynamic_path()
    plot_encoder_hidden_validation()
    plot_sun_fault_diagnosis()
    print(f"Saved core MSSP figures to {OUT}")


if __name__ == "__main__":
    main()
