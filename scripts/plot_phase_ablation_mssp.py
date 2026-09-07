from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np


# Traditional MSSP/Elsevier engineering-paper typography.  Keep SVG text editable.
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 7.0,
        "axes.labelsize": 7.0,
        "axes.titlesize": 7.0,
        "axes.linewidth": 0.55,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "xtick.major.width": 0.50,
        "ytick.major.width": 0.50,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "lines.solid_capstyle": "butt",
    }
)


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = RESULTS / "paper_figures_mssp"
OUT.mkdir(parents=True, exist_ok=True)

METHODS = [
    "raw_nominal",
    "common_phase_only",
    "unconstrained_safe_crossfit",
    "hard3k_all_sensor_fullfit",
    "hard3k_safe_fullfit",
    "hard3k_safe_crossfit",
]

DISPLAY = {
    "raw_nominal": "Raw nominal phase",
    "common_phase_only": "Common phase only",
    "unconstrained_safe_crossfit": "Without three-planet spatial constraint",
    "hard3k_all_sensor_fullfit": "All sensors without topology safeguard",
    "hard3k_safe_fullfit": "Without cross-fitting",
    "hard3k_safe_crossfit": "Full constrained safe method",
}

SPECTRUM_BLUE = "#0000FF"
CYAN = "#00CFEA"
MAGENTA = "#EA00F7"
LIGHT_CYAN = "#77DDE7"
CYCLE_COLORS = [CYAN, MAGENTA, LIGHT_CYAN]
MEAN_COLOR = "#111111"
METHOD_TAG_COLOR = "#D62728"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def number(row: dict[str, str], field: str) -> float:
    return float(row[field])


def rows_by_method(
    rows: Iterable[dict[str, str]], condition: str
) -> dict[str, dict[str, str]]:
    selected = {row["method"]: row for row in rows if row["condition"] == condition}
    missing = set(METHODS) - set(selected)
    if missing:
        raise ValueError(f"Missing methods for {condition}: {sorted(missing)}")
    return selected


def save_bundle(fig: plt.Figure, stem: str) -> list[Path]:
    outputs: list[Path] = []
    for extension in ("svg", "pdf", "png", "tiff"):
        target = OUT / f"{stem}.{extension}"
        kwargs: dict[str, object] = {
            "facecolor": "white",
            "bbox_inches": "tight",
            "pad_inches": 0.025,
        }
        if extension == "png":
            kwargs["dpi"] = 300
        elif extension == "tiff":
            kwargs["dpi"] = 600
        fig.savefig(target, **kwargs)
        outputs.append(target)
    plt.close(fig)
    return outputs


def nearest_value(frequency: np.ndarray, amplitude: np.ndarray, target: float) -> float:
    return float(amplitude[int(np.argmin(np.abs(frequency - target)))])


def make_vertical_spectra(
    metrics: list[dict[str, str]], spectra: list[dict[str, str]]
) -> list[Path]:
    """Six vertically arranged linear envelope spectra in the Xiaochi et al. idiom."""
    fault_metrics = rows_by_method(metrics, "PF50_planet_fault")
    fault_frequency = number(
        fault_metrics["hard3k_safe_crossfit"], "planet_fault_frequency_hz"
    )

    groups: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for method in METHODS:
        selected = [
            row
            for row in spectra
            if row["condition"] == "PF50_planet_fault"
            and row["method"] == method
            and row["signal_component"] == "common"
        ]
        frequency = np.asarray([number(row, "frequency_hz") for row in selected])
        # Milliscale display avoids unreadable 1e-3 tick offsets while preserving linear amplitude.
        amplitude = 1.0e3 * np.asarray([number(row, "amplitude") for row in selected])
        groups[method] = frequency, amplitude

    fig, axes = plt.subplots(
        6,
        1,
        figsize=(130 / 25.4, 238 / 25.4),
        sharex=True,
    )
    fig.subplots_adjust(left=0.135, right=0.985, top=0.985, bottom=0.055, hspace=1.16)

    for index, (ax, method) in enumerate(zip(axes, METHODS)):
        frequency, amplitude = groups[method]
        mask = (frequency >= 0.0) & (frequency <= 60.0)
        frequency_view = frequency[mask]
        amplitude_view = amplitude[mask]
        ax.plot(frequency_view, amplitude_view, color=SPECTRUM_BLUE, lw=0.34)

        y_max = 1.08 * float(np.nanmax(amplitude_view))
        ax.set_xlim(0.0, 60.0)
        ax.set_ylim(0.0, y_max)
        ax.set_xticks(np.arange(0.0, 61.0, 10.0))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
        # Reference MSSP figures retain a complete frequency axis beneath every panel.
        ax.tick_params(direction="out", pad=1.5, labelbottom=True)
        ax.set_ylabel(r"$A\;(\times 10^{-3})$", labelpad=3.0)
        ax.set_xlabel(r"$f\,/\,\mathrm{Hz}$", labelpad=1.0)

        for harmonic in range(1, 5):
            target = harmonic * fault_frequency
            tip = nearest_value(frequency_view, amplitude_view, target)
            text_height = y_max * (0.80 if harmonic % 2 else 0.93)
            harmonic_label = (
                r"$f_{\mathrm{pf}}$"
                if harmonic == 1
                else rf"${harmonic}f_{{\mathrm{{pf}}}}$"
            )
            ax.annotate(
                harmonic_label,
                xy=(target, tip),
                xytext=(target, text_height),
                ha="center",
                va="bottom",
                fontsize=6.4,
                color="black",
                arrowprops={
                    "arrowstyle": "-|>",
                    "mutation_scale": 5.0,
                    "lw": 0.48,
                    "color": "black",
                    "shrinkA": 1.0,
                    "shrinkB": 1.0,
                },
            )

        ax.text(
            0.5,
            -0.50,
            f"({chr(97 + index)})  {DISPLAY[method]}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=6.8,
        )

    return save_bundle(fig, "phase_ablation_envelope_spectra_mssp")


def cycle_matrix(
    rows: list[dict[str, str]], method: str
) -> tuple[np.ndarray, np.ndarray]:
    selected = [row for row in rows if row["method"] == method]
    repeats = sorted({int(float(row["repeat_index"])) for row in selected})
    phase = np.asarray(
        sorted({number(row, "mechanical_repeat_phase_deg") for row in selected})
    )
    phase_index = {value: index for index, value in enumerate(phase)}
    repeat_index = {value: index for index, value in enumerate(repeats)}
    matrix = np.full((len(repeats), len(phase)), np.nan)
    for row in selected:
        matrix[
            repeat_index[int(float(row["repeat_index"]))],
            phase_index[number(row, "mechanical_repeat_phase_deg")],
        ] = number(row, "normalized_amplitude")
    if np.isnan(matrix).any():
        raise ValueError(f"Incomplete cycle matrix for {method}")
    return phase, matrix


def close_curve(theta: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.r_[theta, theta[0] + 2.0 * np.pi], np.r_[values, values[0]]


def rf_graph(
    ax: plt.Axes,
    phase: np.ndarray,
    cycles: np.ndarray,
    radial_scale: float,
    mean_sd: float,
) -> None:
    theta = np.deg2rad(phase)
    for index, trace in enumerate(cycles):
        radius = 1.0 + radial_scale * trace
        theta_closed, radius_closed = close_curve(theta, radius)
        ax.plot(
            theta_closed,
            radius_closed,
            color=CYCLE_COLORS[index % len(CYCLE_COLORS)],
            lw=0.38,
            alpha=0.82,
            zorder=2,
        )

    mean_radius = 1.0 + radial_scale * np.mean(cycles, axis=0)
    theta_closed, mean_closed = close_curve(theta, mean_radius)
    ax.plot(theta_closed, mean_closed, color=MEAN_COLOR, lw=0.62, zorder=3)
    ax.fill_between(np.linspace(0, 2.0 * np.pi, 361), 0.0, 0.64, color="white", zorder=4)
    ax.plot(np.linspace(0, 2.0 * np.pi, 361), np.full(361, 0.64), color="#AA4A44", lw=0.55, zorder=5)
    ax.text(
        0.0,
        0.32,
        f"Mean SD\n{mean_sd:.3f}",
        ha="center",
        va="center",
        fontsize=6.6,
        zorder=6,
    )
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0.0, 1.28)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    ax.spines["polar"].set_visible(False)


def overlay_graph(
    ax: plt.Axes,
    phase: np.ndarray,
    cycles: np.ndarray,
    method_tag: str,
    show_legend: bool,
) -> None:
    for index, trace in enumerate(cycles):
        ax.plot(
            phase,
            trace,
            color=CYCLE_COLORS[index % len(CYCLE_COLORS)],
            lw=0.48,
            alpha=0.72,
            label=f"Repeat {index + 1}",
        )
    ax.plot(
        phase,
        np.mean(cycles, axis=0),
        color=MEAN_COLOR,
        lw=0.80,
        label="Cycle mean",
        zorder=4,
    )
    ax.set_xlim(0.0, 360.0)
    ax.set_xticks(np.arange(0.0, 361.0, 60.0))
    ax.set_xlabel("Mechanical repeat phase (deg)")
    ax.set_ylabel("Normalized amplitude")
    ax.tick_params(direction="out", pad=1.5)
    ax.text(
        0.985,
        0.94,
        method_tag,
        transform=ax.transAxes,
        color=METHOD_TAG_COLOR,
        fontsize=6.6,
        ha="right",
        va="top",
    )
    if show_legend:
        ax.legend(
            loc="upper right",
            bbox_to_anchor=(0.99, 0.80),
            ncol=2,
            fontsize=5.8,
            handlelength=1.5,
            columnspacing=0.8,
            handletextpad=0.35,
        )


def panel_label(ax: plt.Axes, label: str, x: float = -0.06) -> None:
    ax.text(
        x,
        1.02,
        f"({label})",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.2,
        fontstyle="italic",
    )


def make_alignment_plate(cycle_rows: list[dict[str, str]]) -> list[Path]:
    raw_phase, raw_cycles = cycle_matrix(cycle_rows, "raw_nominal")
    full_phase, full_cycles = cycle_matrix(cycle_rows, "hard3k_safe_crossfit")

    raw_sd = float(np.mean(np.std(raw_cycles, axis=0)))
    full_sd = float(np.mean(np.std(full_cycles, axis=0)))
    absolute_max = max(float(np.max(np.abs(raw_cycles))), float(np.max(np.abs(full_cycles))))
    radial_scale = 0.23 / absolute_max

    fig = plt.figure(figsize=(183 / 25.4, 119 / 25.4))
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=[0.78, 1.55],
        left=0.055,
        right=0.985,
        bottom=0.105,
        top=0.965,
        wspace=0.24,
        hspace=0.46,
    )
    ax_a = fig.add_subplot(grid[0, 0], projection="polar")
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0], projection="polar")
    ax_d = fig.add_subplot(grid[1, 1])

    rf_graph(ax_a, raw_phase, raw_cycles, radial_scale, raw_sd)
    overlay_graph(ax_b, raw_phase, raw_cycles, "Raw nominal phase", show_legend=True)
    rf_graph(ax_c, full_phase, full_cycles, radial_scale, full_sd)
    overlay_graph(
        ax_d,
        full_phase,
        full_cycles,
        "Full constrained safe method",
        show_legend=False,
    )

    all_values = np.concatenate([raw_cycles.ravel(), full_cycles.ravel()])
    y_min = float(np.floor((np.min(all_values) - 0.15) * 2.0) / 2.0)
    y_max = float(np.ceil((np.max(all_values) + 0.15) * 2.0) / 2.0)
    for ax in (ax_b, ax_d):
        ax.set_ylim(y_min, y_max)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

    panel_label(ax_a, "a", x=-0.15)
    panel_label(ax_b, "b", x=-0.08)
    panel_label(ax_c, "c", x=-0.15)
    panel_label(ax_d, "d", x=-0.08)

    return save_bundle(fig, "phase_alignment_rf_cycle_mssp")


def compact_path(value: str) -> str:
    return {
        "none": "None",
        "pair_1_2": "Pair 1-2",
        "pair_2_3": "Pair 2-3",
        "complete_graph": "Complete graph",
        "no_path": "No path",
    }.get(value, value)


def write_ablation_table(metrics: list[dict[str, str]]) -> list[Path]:
    healthy = rows_by_method(metrics, "BL600_healthy")
    fault = rows_by_method(metrics, "PF50_planet_fault")
    csv_path = OUT / "phase_ablation_table_mssp.csv"
    fields = [
        "Configuration",
        "Healthy HMPC",
        "PF50 HMPC",
        "PF50 PCC",
        "Fault-family SNR (dB)",
        "Sensors",
        "Selected path",
        "Amplitude error",
    ]
    table_rows: list[dict[str, str]] = []
    for method in METHODS:
        h = healthy[method]
        f = fault[method]
        table_rows.append(
            {
                "Configuration": DISPLAY[method],
                "Healthy HMPC": f"{number(h, 'hmpc'):.3f}",
                "PF50 HMPC": f"{number(f, 'hmpc'):.3f}",
                "PF50 PCC": f"{number(f, 'pcc'):.3f}",
                "Fault-family SNR (dB)": f"{number(f, 'common_fault_family_to_background_db'):.1f}",
                "Sensors": f["sensor_subset"],
                "Selected path": compact_path(f["path_topology"]),
                "Amplitude error": f"{number(f, 'amplitude_preservation_error'):.1e}",
            }
        )

    with csv_path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(table_rows)

    tex_escape = {
        "&": r"\&",
        "%": r"\%",
        "_": r"\_",
    }

    def tex(value: str) -> str:
        for source, target in tex_escape.items():
            value = value.replace(source, target)
        return value

    tex_lines = [
        r"\begin{table*}[t]",
        r"\centering",
        r"\caption{Controlled ablation on the healthy and PF50 measured records.}",
        r"\label{tab:phase_ablation}",
        r"\resizebox{\textwidth}{!}{%",
        r"\begin{tabular}{lccccccc}",
        r"\toprule",
        r"Configuration & Healthy HMPC & PF50 HMPC & PF50 PCC & SNR (dB) & Sensors & Path & Amplitude error \\",
        r"\midrule",
    ]
    for row in table_rows:
        tex_lines.append(
            " & ".join(
                tex(row[field])
                for field in [
                    "Configuration",
                    "Healthy HMPC",
                    "PF50 HMPC",
                    "PF50 PCC",
                    "Fault-family SNR (dB)",
                    "Sensors",
                    "Selected path",
                    "Amplitude error",
                ]
            )
            + r" \\"
        )
    tex_lines.extend(
        [r"\bottomrule", r"\end{tabular}", r"}", r"\end{table*}"]
    )
    tex_path = OUT / "phase_ablation_table_mssp.tex"
    tex_path.write_text("\n".join(tex_lines) + "\n", encoding="utf-8")
    return [csv_path, tex_path]


def write_notes() -> list[Path]:
    legend = """# Draft MSSP-style figure legends

## Six-configuration common-mode envelope spectra
Common-mode envelope spectra of the measured PF50 record under six controlled configurations. Black arrows indicate the first four theoretical planet-fault harmonics, where $f_{pf}=10.757$ Hz. Linear spectral amplitude is displayed in milliscale units. The panels use independent vertical limits, following conventional envelope-spectrum presentation; quantitative comparisons are reported in the accompanying ablation table.

## Phase-alignment RF and cycle comparison
Cycle consistency before and after the proposed phase correction on the measured PF50 record. (a,c) RF-style circular views of the same three complete 31-carrier-turn mechanical repeats. Radial deviation represents normalized vibration amplitude with one common scale for both panels and is not a physical orbit. The centre reports the phase-wise mean cycle standard deviation. (b,d) The corresponding cycle overlays; coloured curves denote individual repeats and the black curve denotes their mean. No error bars or significance test are shown because one 50 s record is available for each condition.
"""
    qa = """# MSSP-style figure QA

- Backend: Python/Matplotlib only for plotting, previewing, and all exports.
- Typography: Times New Roman serif family; STIX fallback and STIX math.
- Source data: exported phase-ablation CSV files; no synthetic data.
- Statistics: one 50 s record per condition and three complete mechanical repeats; no fabricated error bars or p-values.
- Spectrum: linear amplitude, independent panel limits stated explicitly, no smoothing or selective peak enhancement.
- RF graph: a normalized circular visualization, not a shaft orbit; one common radial scale is used before and after correction.
- SVG text remains editable; PDF uses TrueType fonts; TIFF is exported at 600 dpi.
"""
    legend_path = OUT / "figure_legends_mssp.md"
    qa_path = OUT / "qa_notes_mssp.md"
    legend_path.write_text(legend, encoding="utf-8")
    qa_path.write_text(qa, encoding="utf-8")
    return [legend_path, qa_path]


def validate_outputs(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"Missing or empty outputs: {missing}")
    for path in paths:
        if path.suffix.lower() == ".svg":
            svg = path.read_text(encoding="utf-8")
            if "<text" not in svg:
                raise RuntimeError(f"SVG text was converted to paths: {path}")


def main() -> None:
    metrics = read_rows(RESULTS / "phase_ablation_metrics.csv")
    spectra = read_rows(RESULTS / "phase_ablation_spectra_source.csv")
    cycles = read_rows(RESULTS / "phase_ablation_cycle_source.csv")

    outputs: list[Path] = []
    outputs.extend(make_vertical_spectra(metrics, spectra))
    outputs.extend(make_alignment_plate(cycles))
    outputs.extend(write_ablation_table(metrics))
    outputs.extend(write_notes())
    validate_outputs(outputs)
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
