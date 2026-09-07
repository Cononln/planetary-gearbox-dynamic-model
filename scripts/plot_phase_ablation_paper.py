from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# Mandatory publication settings: sans-serif font and editable SVG text.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7.0
plt.rcParams["axes.linewidth"] = 0.75
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["legend.frameon"] = False
plt.rcParams["xtick.major.width"] = 0.65
plt.rcParams["ytick.major.width"] = 0.65
plt.rcParams["xtick.major.size"] = 3.0
plt.rcParams["ytick.major.size"] = 3.0


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = RESULTS / "paper_figures"
OUT.mkdir(parents=True, exist_ok=True)

METHODS = [
    "raw_nominal",
    "common_phase_only",
    "unconstrained_safe_crossfit",
    "hard3k_all_sensor_fullfit",
    "hard3k_safe_fullfit",
    "hard3k_safe_crossfit",
]
SHORT = {
    "raw_nominal": "Raw",
    "common_phase_only": "Common",
    "unconstrained_safe_crossfit": "No 3k",
    "hard3k_all_sensor_fullfit": "All 3 ch.",
    "hard3k_safe_fullfit": "No CF",
    "hard3k_safe_crossfit": "Full",
}
LONG = {
    "raw_nominal": "Raw nominal phase",
    "common_phase_only": "Common phase only",
    "unconstrained_safe_crossfit": "Without three-planet spatial constraint",
    "hard3k_all_sensor_fullfit": "All sensors without topology safeguard",
    "hard3k_safe_fullfit": "Without cross-fitting",
    "hard3k_safe_crossfit": "Full constrained safe method",
}
COLORS = {
    "raw_nominal": "#D5D5D5",
    "common_phase_only": "#B9CEE5",
    "unconstrained_safe_crossfit": "#88AAD0",
    "hard3k_all_sensor_fullfit": "#5F8DBB",
    "hard3k_safe_fullfit": "#376F9F",
    "hard3k_safe_crossfit": "#0F4D92",
}
HATCHES = ["//", "\\\\", "..", "xx", "++", ""]
HEALTHY_COLOR = "#767676"
FAULT_COLOR = "#B64342"
FAULT_LINE_COLOR = "#B64342"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def number(row: dict[str, str], name: str) -> float:
    return float(row[name])


def rows_by_method(
    rows: Iterable[dict[str, str]], condition: str
) -> dict[str, dict[str, str]]:
    result = {}
    for row in rows:
        if row["condition"] == condition:
            result[row["method"]] = row
    missing = set(METHODS) - set(result)
    if missing:
        raise ValueError(f"Missing methods for {condition}: {sorted(missing)}")
    return result


def style_axis(ax: plt.Axes) -> None:
    ax.tick_params(direction="out", pad=2)
    ax.grid(False)


def panel_label(ax: plt.Axes, label: str, x: float = -0.13) -> None:
    ax.text(
        x,
        1.035,
        label,
        transform=ax.transAxes,
        fontsize=8.5,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def annotate_vertical_bars(
    ax: plt.Axes, bars, values: list[float], decimals: int = 2
) -> None:
    lo, hi = ax.get_ylim()
    offset = 0.018 * (hi - lo)
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + offset,
            f"{value:.{decimals}f}",
            ha="center",
            va="bottom",
            fontsize=5.8,
        )


def method_bars(
    ax: plt.Axes,
    rows: dict[str, dict[str, str]],
    metric: str,
    ylabel: str,
    ylim: tuple[float, float],
    decimals: int,
) -> None:
    values = [number(rows[m], metric) for m in METHODS]
    x = np.arange(len(METHODS))
    bars = ax.bar(
        x,
        values,
        width=0.72,
        color=[COLORS[m] for m in METHODS],
        edgecolor="#272727",
        linewidth=0.55,
    )
    for bar, hatch in zip(bars, HATCHES):
        bar.set_hatch(hatch)
    ax.set_xticks(x, [SHORT[m] for m in METHODS], rotation=24, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_ylim(*ylim)
    annotate_vertical_bars(ax, bars, values, decimals)
    style_axis(ax)


def save_bundle(fig: plt.Figure, stem: str) -> list[Path]:
    outputs = []
    for extension in ("svg", "pdf", "png", "tiff"):
        path = OUT / f"{stem}.{extension}"
        kwargs = {}
        if extension == "png":
            kwargs["dpi"] = 300
        elif extension == "tiff":
            kwargs["dpi"] = 600
        fig.savefig(path, facecolor="white", **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def make_summary(metrics: list[dict[str, str]]) -> list[Path]:
    healthy = rows_by_method(metrics, "BL600_healthy")
    fault = rows_by_method(metrics, "PF50_planet_fault")

    fig = plt.figure(figsize=(183 / 25.4, 128 / 25.4), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.08])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    method_bars(ax_a, fault, "hmpc", "HMPC", (0.0, 0.56), 3)
    ax_a.set_title("Held-out phase concentration", pad=4)
    panel_label(ax_a, "a")

    method_bars(ax_b, fault, "pcc", "PCC", (0.0, 0.78), 3)
    ax_b.set_title("Waveform consistency", pad=4)
    panel_label(ax_b, "b")

    method_bars(
        ax_c,
        fault,
        "common_fault_family_to_background_db",
        "Common-mode fault family (dB)",
        (18.0, 26.0),
        1,
    )
    ax_c.set_title("Planet-fault feature visibility", pad=4)
    panel_label(ax_c, "c")

    y = np.arange(len(METHODS))[::-1]
    healthy_values = np.array([number(healthy[m], "hmpc") for m in METHODS])
    fault_values = np.array([number(fault[m], "hmpc") for m in METHODS])
    for yi, h, f in zip(y, healthy_values, fault_values):
        ax_d.plot([h, f], [yi, yi], color="#BEBEBE", lw=0.8, zorder=1)
    ax_d.scatter(
        healthy_values,
        y,
        s=22,
        facecolors="white",
        edgecolors=HEALTHY_COLOR,
        linewidths=0.8,
        label="BL600 healthy",
        zorder=3,
    )
    ax_d.scatter(
        fault_values,
        y,
        s=24,
        color=FAULT_COLOR,
        edgecolors="white",
        linewidths=0.4,
        label="PF50 planet fault",
        zorder=4,
    )
    ax_d.set_yticks(y, [SHORT[m] for m in METHODS])
    ax_d.set_xlabel("HMPC")
    ax_d.set_xlim(0.05, 0.64)
    ax_d.set_title("Healthy/fault response of each ablation", pad=4)
    ax_d.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        ncol=2,
        fontsize=5.8,
        handletextpad=0.45,
        columnspacing=1.2,
    )
    full_y = y[METHODS.index("hard3k_safe_crossfit")]
    ax_d.text(
        healthy_values[-1] - 0.008,
        full_y + 0.16,
        "no path",
        color=HEALTHY_COLOR,
        fontsize=5.3,
        ha="right",
    )
    ax_d.text(
        fault_values[-1] + 0.008,
        full_y + 0.16,
        "pair 1-2",
        color=FAULT_COLOR,
        fontsize=5.3,
        ha="left",
    )
    style_axis(ax_d)
    panel_label(ax_d, "d")

    return save_bundle(fig, "phase_ablation_summary")


def make_spectra(
    metrics: list[dict[str, str]], spectra: list[dict[str, str]]
) -> list[Path]:
    fault_metrics = rows_by_method(metrics, "PF50_planet_fault")
    groups: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for method in METHODS:
        selected = [
            row
            for row in spectra
            if row["condition"] == "PF50_planet_fault"
            and row["method"] == method
            and row["signal_component"] == "common"
        ]
        groups[method] = (
            np.array([number(row, "frequency_hz") for row in selected]),
            np.array([number(row, "relative_to_background_db") for row in selected]),
        )

    all_db = np.concatenate([groups[m][1] for m in METHODS])
    y_low = max(-20.0, float(np.floor(np.nanpercentile(all_db, 1) / 5) * 5))
    y_high = min(70.0, float(np.ceil(np.nanpercentile(all_db, 99.8) / 5) * 5))
    fault_frequency = number(fault_metrics["hard3k_safe_crossfit"], "planet_fault_frequency_hz")

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(183 / 25.4, 116 / 25.4),
        sharex=True,
        sharey=True,
        layout="constrained",
    )
    for index, (ax, method) in enumerate(zip(axes.flat, METHODS)):
        frequency, amplitude_db = groups[method]
        ax.plot(frequency, amplitude_db, color="#235A91", lw=0.55)
        for harmonic in range(1, 5):
            target = harmonic * fault_frequency
            ax.axvline(
                target,
                color=FAULT_LINE_COLOR,
                lw=0.55,
                ls=(0, (2, 2)),
                alpha=0.85,
            )
            if index == 0:
                ax.text(
                    target,
                    y_high - 1.5,
                    rf"${harmonic}f_{{pf}}$" if harmonic > 1 else r"$f_{pf}$",
                    color=FAULT_LINE_COLOR,
                    fontsize=5.6,
                    ha="center",
                    va="top",
                )
        snr = number(
            fault_metrics[method], "common_fault_family_to_background_db"
        )
        ax.text(
            0.97,
            0.94,
            f"SNR = {snr:.1f} dB",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=5.8,
            bbox={"boxstyle": "round,pad=0.18", "fc": "white", "ec": "#BEBEBE", "lw": 0.45},
        )
        ax.set_title(f"({chr(97 + index)})  {LONG[method]}", loc="left", pad=4)
        ax.set_xlim(0, 60)
        ax.set_ylim(y_low, y_high)
        style_axis(ax)

    fig.supxlabel("Envelope frequency (Hz)", fontsize=7.2)
    fig.supylabel("Common-mode amplitude / background (dB)", fontsize=7.2)
    return save_bundle(fig, "phase_ablation_common_mode_spectra")


def matrix_from_rows(
    rows: list[dict[str, str]],
    method: str,
    row_field: str,
    column_field: str,
    value_field: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    selected = [row for row in rows if row["method"] == method]
    row_values = np.array(sorted({number(row, row_field) for row in selected}))
    column_values = np.array(sorted({number(row, column_field) for row in selected}))
    row_index = {value: i for i, value in enumerate(row_values)}
    column_index = {value: i for i, value in enumerate(column_values)}
    matrix = np.full((len(row_values), len(column_values)), np.nan)
    for row in selected:
        matrix[row_index[number(row, row_field)], column_index[number(row, column_field)]] = number(
            row, value_field
        )
    return row_values, column_values, matrix


def make_consistency(
    map_rows: list[dict[str, str]], cycle_rows: list[dict[str, str]]
) -> list[Path]:
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(183 / 25.4, 110 / 25.4),
        layout="constrained",
    )
    methods = ["raw_nominal", "hard3k_safe_crossfit"]
    titles = ["Raw nominal phase", "Full constrained method"]
    map_images = []
    cycle_matrices = []

    for column, (method, title) in enumerate(zip(methods, titles)):
        harmonic, phase, matrix = matrix_from_rows(
            map_rows,
            method,
            "mesh_harmonic",
            "mechanical_repeat_phase_deg",
            "phase_concentration",
        )
        image = axes[0, column].imshow(
            matrix,
            origin="lower",
            aspect="auto",
            extent=[phase.min(), phase.max(), harmonic.min() - 0.5, harmonic.max() + 0.5],
            cmap="Blues",
            vmin=0,
            vmax=1,
            interpolation="nearest",
        )
        map_images.append(image)
        axes[0, column].set_title(title, pad=4)
        axes[0, column].set_xlabel("Mechanical repeat phase (deg)")
        axes[0, column].set_yticks([1, 2, 3, 4])
        if column == 0:
            axes[0, column].set_ylabel("Mesh harmonic")
        style_axis(axes[0, column])
        panel_label(axes[0, column], chr(97 + column), x=-0.10)

        repeat, cycle_phase, cycle_matrix = matrix_from_rows(
            cycle_rows,
            method,
            "repeat_index",
            "mechanical_repeat_phase_deg",
            "normalized_amplitude",
        )
        cycle_matrices.append((cycle_phase, cycle_matrix))
        mean_cycle = np.mean(cycle_matrix, axis=0)
        std_cycle = np.std(cycle_matrix, axis=0)
        for repeat_index, trace in enumerate(cycle_matrix, start=1):
            axes[1, column].plot(
                cycle_phase,
                trace,
                color="#8FAFCF",
                lw=0.35,
                alpha=0.42,
                label="Individual repeats" if repeat_index == 1 else None,
            )
        axes[1, column].fill_between(
            cycle_phase,
            mean_cycle - std_cycle,
            mean_cycle + std_cycle,
            color="#B9CEE5",
            alpha=0.55,
            linewidth=0,
            label="Mean +/- 1 SD",
        )
        axes[1, column].plot(
            cycle_phase,
            mean_cycle,
            color="#0F4D92",
            lw=0.75,
            label="Cycle mean",
        )
        axes[1, column].set_title(f"Cycle overlay, {title.lower()}", pad=4)
        axes[1, column].set_xlabel("Mechanical repeat phase (deg)")
        if column == 0:
            axes[1, column].set_ylabel("Normalized amplitude")
        axes[1, column].text(
            0.98,
            0.94,
            f"mean SD = {np.mean(std_cycle):.3f}",
            transform=axes[1, column].transAxes,
            ha="right",
            va="top",
            fontsize=5.6,
            bbox={"boxstyle": "round,pad=0.16", "fc": "white", "ec": "#BEBEBE", "lw": 0.4},
        )
        style_axis(axes[1, column])
        panel_label(axes[1, column], chr(99 + column), x=-0.10)

    colorbar_map = fig.colorbar(map_images[-1], ax=axes[0, :], shrink=0.85, pad=0.025)
    colorbar_map.set_label("HMPC phase concentration")
    y_min = min(float(np.nanpercentile(matrix, 0.5)) for _, matrix in cycle_matrices)
    y_max = max(float(np.nanpercentile(matrix, 99.5)) for _, matrix in cycle_matrices)
    margin = 0.08 * (y_max - y_min)
    for ax in axes[1, :]:
        ax.set_ylim(y_min - margin, y_max + margin)
    axes[1, 1].legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        ncol=3,
        fontsize=5.5,
        handlelength=1.6,
    )
    return save_bundle(fig, "phase_alignment_consistency")


def validate_outputs(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"Missing or empty figure outputs: {missing}")
    for path in paths:
        if path.suffix == ".svg":
            text = path.read_text(encoding="utf-8")
            if "<text" not in text:
                raise RuntimeError(f"SVG text is not editable: {path}")


def write_documentation() -> None:
    contract = """# Figure contract

Core conclusion: The full method separates common speed phase from moving-path phase while the three-planet constraint and held-out topology safeguard prevent unsupported correction of healthy data.

Figure archetype: quantitative grid plus before/after phase-consistency plate.

Target/output: double-column high-impact signal-processing journal figure; editable SVG/PDF, 600 dpi TIFF, and PNG preview.

Backend: Python/Matplotlib only.

Final size: 183 mm wide; 110-128 mm high.

Panel evidence:
- Ablation summary: held-out HMPC, PCC, fault-family visibility, and healthy/fault response.
- Common-mode spectra: first four theoretical planet-fault harmonics for all six configurations.
- Consistency plate: raw versus full HMPC maps and three complete mechanical-repeat cycles.

Statistics: one 50 s record per condition; three complete 31-carrier-turn repeats for HMPC/PCC; no independent-record significance test or fabricated error bars.

Reviewer risk: the record count limits population-level inference; the figures support controlled within-record ablation and mechanism validation only.
"""
    legends = """# Draft figure legends

## Phase-ablation summary
Controlled ablation on the measured PF50 record. (a) Held-out multi-channel phase concentration (HMPC). (b) Phase consistency coefficient (PCC). (c) Signal-to-background ratio of the first four theoretical planet-fault lines in the aligned common-mode envelope. (d) Healthy and fault HMPC for each configuration. Bars show deterministic values from one 50 s record; HMPC and PCC use three complete 31-carrier-turn mechanical repeats. CF, cross-fitting; 3k, three-planet spatial-order constraint.

## Common-mode envelope spectra
Common-mode envelope spectra of the PF50 record for six controlled configurations. Red dashed lines mark the first four theoretical planet-fault harmonics. Each spectrum is expressed relative to its own off-line 2-60 Hz median background; the inset reports the corresponding four-line RMS signal-to-background ratio.

## Phase-alignment consistency
Raw nominal phase versus the full constrained safe method for the PF50 record. (a,b) HMPC as a function of mechanical-repeat phase and mesh harmonic. (c,d) Normalized synchronous cycle matrices for the three available complete 31-carrier-turn repeats. The same color limits are used within each row.
"""
    qa = """# Figure QA notes

- All panels were generated by Python/Matplotlib from exported CSV source data.
- SVG text remains editable and PDF fonts use TrueType embedding.
- PNG previews are 300 dpi; TIFF exports are 600 dpi.
- Method colors are consistent across figures and hatch patterns support grayscale printing.
- Axes that invite direct comparison share limits.
- No error bars or p-values are shown because only one record per condition is available.
- No raster image enhancement, selective smoothing, or local contrast manipulation was applied.
"""
    (OUT / "figure_contract.md").write_text(contract, encoding="utf-8")
    (OUT / "figure_legends.md").write_text(legends, encoding="utf-8")
    (OUT / "qa_notes.md").write_text(qa, encoding="utf-8")


def main() -> None:
    metrics = read_rows(RESULTS / "phase_ablation_metrics.csv")
    spectra = read_rows(RESULTS / "phase_ablation_spectra_source.csv")
    map_rows = read_rows(RESULTS / "phase_ablation_hmpc_map_source.csv")
    cycle_rows = read_rows(RESULTS / "phase_ablation_cycle_source.csv")

    outputs = []
    outputs.extend(make_summary(metrics))
    outputs.extend(make_spectra(metrics, spectra))
    outputs.extend(make_consistency(map_rows, cycle_rows))
    validate_outputs(outputs)
    write_documentation()
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
