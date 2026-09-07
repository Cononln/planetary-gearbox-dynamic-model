"""Measured feasibility figures: first-pass evidence and exploratory diagnosis.

All panels use exported real measurements. Dots are 10 s blocks from ONE
record per condition; no independent-replicate significance is inferred.
"""
from __future__ import annotations

import csv
import os
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "runtime" / "phase_python"))
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "runtime" / "matplotlib_feasibility"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = ROOT / "results" / "healthy_reference_feasibility_20260905"
OUT = BASE / "figures"
OUT.mkdir(exist_ok=True)
RUNS = {
    "initial": BASE,
    "anchored": Path(str(BASE) + "_exploratory_calibration_anchored"),
    "nominal": Path(str(BASE) + "_exploratory_nominal"),
}
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": .7, "lines.linewidth": 1.1,
    "xtick.major.width": .7, "ytick.major.width": .7,
    "legend.frameon": False, "legend.fontsize": 6,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
COL = {"BL": "#8195A9", "PF50": "#B56E55", "raw": "#6C7885", "frozen": "#355F83", "diff": "#AF8A50"}


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def values(rows, field, **filters):
    return np.array([float(r[field]) for r in rows if all(r[k] == str(v) for k, v in filters.items())])


DATA = {key: read(path / "heldout_metrics.csv") for key, path in RUNS.items()}
ENC = {key: read(path / "encoder_validation.csv") for key, path in RUNS.items()}
TRACE = {key: read(path / "angle_trace.csv") for key, path in RUNS.items()}


def mark(ax, label, title):
    ax.text(-.17, 1.06, label, transform=ax.transAxes, fontweight="bold", fontsize=9)
    ax.set_title(title, loc="left", pad=8)


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=220)
    fig.savefig(OUT / f"{name}.svg")
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    svg = ET.parse(OUT / f"{name}.svg")
    assert len(svg.findall(".//{http://www.w3.org/2000/svg}text")) > 15
    plt.close(fig)


def summarize():
    rows = []
    methods = list(dict.fromkeys(r["method"] for r in DATA["initial"]))
    for run, data in DATA.items():
        for condition in ("BL", "PF50"):
            for method in methods:
                row = {"run": run, "condition": condition, "method": method, "n_time_blocks": 3}
                for metric in ("coherence", "fp_family_db", "twofp_family_db", "retained_diff_fp_db", "retained_diff_twofp_db", "amplitude_error"):
                    v = values(data, metric, condition=condition, method=method)
                    assert len(v) == 3 and np.isfinite(v).all()
                    row[metric + "_mean"] = v.mean()
                    row[metric + "_min"] = v.min()
                    row[metric + "_max"] = v.max()
                rows.append(row)
    with (BASE / "summary.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    for run in RUNS:
        print("\nRUN", run)
        for condition in ("BL", "PF50"):
            for m in ("common_all", "common_selected", "healthy_constant", "healthy_frozen_all", "healthy_frozen_selected", "within_record", "sensor1"):
                row = next(r for r in rows if (r["run"], r["condition"], r["method"]) == (run, condition, m))
                print(condition, m, "C %.6f fp %.4f 2fp %.4f" % tuple(row[f + "_mean"] for f in ("coherence", "fp_family_db", "twofp_family_db")))
            print(condition, "angle_RMSE_blocks", values(ENC[run], "carrier_angle_rmse_deg", condition=condition))
    return rows


def first_pass():
    fig, axes = plt.subplots(2, 2, figsize=(183 / 25.4, 132 / 25.4))
    fig.subplots_adjust(left=.11, right=.98, bottom=.14, top=.90, wspace=.40, hspace=.62)
    fig.suptitle("Frozen healthy-reference phase correction: first-pass held-out results", fontsize=10, y=.985)
    ax = axes[0, 0]
    for condition in ("BL", "PF50"):
        for block in (1, 2, 3):
            x = values(TRACE["initial"], "time_s", condition=condition, block=block)
            y = values(TRACE["initial"], "estimated_error_deg", condition=condition, block=block)
            ax.plot(x, y, color=COL[condition], label=condition if block == 1 else None)
    ax.axhline(0, color=".7", lw=.6, ls="--")
    ax.set(xlabel="Time (s)", ylabel="Carrier-angle error (deg)", xlim=(27, 59))
    ax.legend(loc="upper left"); mark(ax, "a", "Common-phase drift vs encoder")
    ax = axes[0, 1]
    methods = ["common_selected", "healthy_constant", "healthy_frozen_selected", "within_record"]
    labels = ["Common\nonly", "Constant\nhealthy", "Periodic\nhealthy", "Within-\nrecord"]
    for condition, offset, marker in [("BL", -.10, "o"), ("PF50", .10, "s")]:
        for i, method in enumerate(methods):
            v = values(DATA["initial"], "coherence", condition=condition, method=method)
            x = i + offset
            ax.scatter(x + np.linspace(-.032, .032, 3), v, s=13, marker=marker, color=COL[condition], alpha=.8)
            ax.plot([x - .075, x + .075], [v.mean()] * 2, color=COL[condition], lw=1.8)
    ax.set(xticks=range(4), xticklabels=labels, ylim=(0, 1.03), ylabel="Inter-channel phase coherence")
    mark(ax, "b", "Same sensor pair (channels 2 + 3)")
    ax = axes[1, 0]
    fields = ["fp_family_db", "twofp_family_db"]
    for metric, offset, marker, color, label in [(fields[0], -.10, "o", "#355F83", "$f_p$ family"), (fields[1], .10, "s", "#AF8A50", "$2f_p$ family")]:
        for i, method in enumerate(methods):
            v = values(DATA["initial"], metric, condition="PF50", method=method)
            ax.scatter(i + offset + np.linspace(-.032, .032, 3), v, s=13, marker=marker, color=color, label=label if i == 0 else None)
            ax.plot([i + offset - .075, i + offset + .075], [v.mean()] * 2, color=color, lw=1.8)
    ax.set(xticks=range(4), xticklabels=labels, ylim=(0, 23), ylabel="Fault-family / background (dB)")
    ax.legend(loc="upper right"); mark(ax, "c", "PF50: coherence is not diagnosis")
    ax = axes[1, 1]
    spectra = read(BASE / "envelope_spectra.csv")
    for method, color, label, style in [("common_selected", COL["raw"], "Common only (2 + 3)", "-"), ("healthy_frozen_selected", COL["frozen"], "Frozen healthy (2 + 3)", "-"), ("retained_raw_differential", COL["diff"], "Retained raw differential (all 3)", "--")]:
        x = values(spectra, "carrier_order", condition="PF50", block=1, method=method)
        y = values(spectra, "amplitude", condition="PF50", block=1, method=method)
        use = (x >= 1) & (x <= 12)
        ax.plot(x[use], y[use], color=color, lw=.8, label=label, ls=style)
    for i in range(1, 5):
        ax.axvline(i * 84 / 31, color=".6", lw=.5, ls=":", zorder=0)
    ax.set(xlabel="Carrier order", ylabel="Normalized envelope amplitude", xlim=(1, 12), ylim=(0, None))
    ax.legend(loc="upper right", fontsize=5.7); mark(ax, "d", "PF50 envelope spectrum, 27-37 s")
    fig.text(.11, .027, "Dots: three 10 s blocks from one record per condition; horizontal ticks: means. No independent-repeat significance test.\nEncoder is used only for validation; one constant angular offset is removed. Spectra share healthy-derived scaling.", fontsize=6, color=".3")
    save(fig, "first_pass_validation")


def exploratory():
    fig, axes = plt.subplots(2, 2, figsize=(183 / 25.4, 132 / 25.4))
    fig.subplots_adjust(left=.11, right=.98, bottom=.14, top=.90, wspace=.40, hspace=.62)
    fig.suptitle("Exploratory failure isolation: frequency anchoring and retained phase residuals", fontsize=9.5, y=.985)
    ax = axes[0, 0]
    for run, color, label, style in [("anchored", "#355F83", "Calibration anchored", "-"), ("nominal", "#AF8A50", "Constant-speed phase", "--")]:
        for block in (1, 2, 3):
            x = values(TRACE[run], "time_s", condition="PF50", block=block)
            y = values(TRACE[run], "estimated_error_deg", condition="PF50", block=block)
            ax.plot(x, y, color=color, lw=.8, label=label if block == 1 else None, ls=style)
    ax.axhline(0, color=".7", lw=.6, ls=":")
    ax.set(xlabel="Time (s)", ylabel="Carrier-angle error (deg)", xlim=(27, 59))
    ax.legend(loc="best"); mark(ax, "a", "PF50: reduced cumulative drift")
    ax = axes[0, 1]
    methods = ["healthy_constant", "healthy_frozen_selected", "within_record"]
    for condition, offset, marker in [("BL", -.08, "o"), ("PF50", .08, "s")]:
        for i, method in enumerate(methods):
            v = values(DATA["anchored"], "coherence", condition=condition, method=method)
            ax.scatter(i + offset + np.linspace(-.025, .025, 3), v, color=COL[condition], marker=marker, s=16, label=condition if i == 0 else None)
            ax.plot([i + offset - .07, i + offset + .07], [v.mean()] * 2, color=COL[condition], lw=1.8)
    ax.set(xticks=range(3), xticklabels=["Constant\nhealthy", "Periodic\nhealthy", "Within-\nrecord"], ylim=(.80, 1), ylabel="Inter-channel phase coherence")
    ax.legend(loc="lower right"); mark(ax, "b", "Periodic gain transfers incompletely")
    ax = axes[1, 0]
    methods = ["common_selected", "healthy_frozen_selected", "within_record", "sensor1"]
    for i, method in enumerate(methods):
        v = values(DATA["anchored"], "fp_family_db", condition="PF50", method=method)
        color = COL["frozen"] if i in (1, 2) else COL["raw"]
        ax.scatter(i + np.linspace(-.04, .04, 3), v, s=16, color=color)
        ax.plot([i - .09, i + .09], [v.mean()] * 2, lw=1.8, color=color)
    ax.set(xticks=range(4), xticklabels=["Common\n(2 + 3)", "Healthy\nfrozen", "Within-\nrecord", "Single\nchannel 1"], ylim=(0, 23), ylabel="$f_p$ family / background (dB)")
    mark(ax, "c", "PF50: feature loss persists")
    ax = axes[1, 1]
    for i, condition in enumerate(("BL", "PF50")):
        v = 1 - values(DATA["anchored"], "coherence", condition=condition, method="healthy_frozen_selected")
        ax.scatter(i + np.linspace(-.04, .04, 3), v, color=COL[condition], s=20, marker="o" if i == 0 else "s")
        ax.plot([i - .1, i + .1], [v.mean()] * 2, color=COL[condition], lw=1.8)
    ax.set(xticks=[0, 1], xticklabels=["Healthy", "PF50"], xlim=(-.5, 1.5), ylim=(0, .15), ylabel="Phase residual (1 - C)")
    mark(ax, "d", "Candidate residual, not classifier proof")
    fig.text(.11, .027, "Exploratory: these test windows had already been inspected. Calibration uses only 3-15 s vibration data; no encoder input.\nDots are within-record blocks, not independent specimens. A phase residual may also reflect transfer-path or mounting changes.", fontsize=6, color=".3")
    save(fig, "exploratory_failure_isolation")


if __name__ == "__main__":
    summarize()
    first_pass()
    exploratory()
    print("\nFigures exported; SVG text nodes verified. Output:", OUT)
