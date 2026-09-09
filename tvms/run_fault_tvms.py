"""Parameterized local root-crack TVMS, without a dynamics simulation.

The straight root-crack geometry follows Peng Yue (2025): a crack starts at
the tooth root, has extension ``q_c`` and forms angle ``gamma`` with the tooth
centreline.  Only the named fault tooth is modified when it is an active pair;
every other active pair is evaluated with the unchanged healthy potential-
energy model.

This module deliberately does *not* import any dynamic response runner.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

# Must precede every Matplotlib import: the global Windows font cache can be
# locked by another active plotting process.
ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "results" / "tvms" / ".matplotlib_fault"))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from scipy.integrate import quad
from scipy.optimize import brentq

import run_tvms as healthy


TVMS_OUT = ROOT / "results" / "tvms"
CONFIG = ROOT / "configs" / "fault_tvms_config.json"
TIMES = FontProperties(fname=r"C:\Windows\Fonts\times.ttf", size=9.5)
TIMES_BOLD = FontProperties(fname=r"C:\Windows\Fonts\timesbd.ttf", size=10.5)
SIMSUN = FontProperties(fname=r"C:\Windows\Fonts\simsun.ttc", size=9.5)


@dataclass(frozen=True)
class Context:
    sun: healthy.Gear
    planet: healthy.Gear
    ring: healthy.Gear
    sp_path: healthy.MeshPath
    pr_path: healthy.MeshPath
    alpha_nom: float
    alpha_sp: float
    alpha_pr: float
    face_sp: float
    face_pr: float
    young: float
    poisson: float
    x_each: float
    ring_foundation_compliance: float
    sp_origin: float
    pr_origin: float
    zs: int
    zp: int
    zr: int
    n_planet: int


def read_config() -> dict:
    return json.loads(CONFIG.read_text(encoding="utf-8"))


def build_context() -> Context:
    value = healthy.source_value
    zs, zp, zr = (int(value(item)) for item in ("gear.zs", "gear.zp", "gear.zr"))
    module, alpha_nom, centre = (value(item) for item in
                                 ("gear.module", "gear.alpha_nominal", "gear.carrierRadius"))
    face_sp, face_pr = (value(item) for item in ("gear.faceWidthSunPlanet", "gear.faceWidthRing"))
    young, poisson = (value(item) for item in ("material.gear.E", "material.gear.nu"))
    sun = healthy.Gear("sun", zs, module * zs / 2, module * zs * math.cos(alpha_nom) / 2,
                       value("gear.addendumRadiusSun"), value("gear.rootRadiusSun"))
    planet = healthy.Gear("planet", zp, module * zp / 2, module * zp * math.cos(alpha_nom) / 2,
                          value("gear.addendumRadiusPlanet"), value("gear.rootRadiusPlanet"))
    ring = healthy.Gear("ring", zr, module * zr / 2, module * zr * math.cos(alpha_nom) / 2,
                        value("gear.addendumRadiusRing"), value("gear.rootRadiusRing"), True)
    alpha_sp = math.acos((sun.pitch + planet.pitch) * math.cos(alpha_nom) / centre)
    alpha_pr = math.acos((ring.pitch - planet.pitch) * math.cos(alpha_nom) / centre)
    base_pitch = math.pi * module * math.cos(alpha_nom)
    sp_path = healthy.external_path(sun, planet, alpha_sp, base_pitch)
    pr_path = healthy.internal_path(planet, ring, alpha_pr, centre, base_pitch)
    x_total = (zs + zp) * (healthy.involute(alpha_sp) - healthy.involute(alpha_nom)) / (2 * math.tan(alpha_nom))

    # This is the frozen healthy PR foundation term.  It is recomputed from
    # the same FEM geometry and mesh resolution solely to make this runner
    # self-contained; it is not changed by q or by tooth identity.
    root_thickness = healthy.section_width(ring, ring.root, alpha_nom)
    result = healthy.annular_rim_normal_compliance(
        root_radius=ring.root, outer_radius=value("gear.outerRadiusRing"),
        face_width=face_pr, young=young, poisson=poisson,
        hole_radius=value("gear.ringHoleRadius"), hole_diameter=value("gear.ringHoleDiameter"),
        hole_count=int(value("gear.ringHoleCount")), tooth_count=zr,
        root_tooth_thickness=root_thickness, pressure_angle=alpha_pr,
        nr=36, ntheta=672,
    )
    return Context(sun, planet, ring, sp_path, pr_path, alpha_nom, alpha_sp, alpha_pr,
                   face_sp, face_pr, young, poisson, x_total / 2, result.compliance_mean,
                   sp_path.epsilon / 2, pr_path.epsilon / 2, zs, zp, zr,
                   int(value("model.nPlanet")))


def root_crack_tooth_compliance(
        gear: healthy.Gear, r_contact: float, face_width: float, alpha_work: float,
        context: Context, crack_depth: float, crack_angle: float, x_equiv: float,
        foundation_compliance: float | None, min_ligament_fraction: float) -> float:
    """Unit-load compliance of one cracked tooth using Peng Yue's section law.

    Let ``h_x`` be half of the healthy local section thickness, ``h_c`` the
    root-section half thickness, and ``h_q=h_c-q_c sin(gamma)``.  Per Eqs.
    (3-1)--(3-2), the remaining section is ``2 h_x`` for ``h_x <= h_q`` and
    ``h_c+h_x-q_c sin(gamma)`` otherwise.  Consequently only bending and
    shear use the cracked I_x/A_x; axial and foundation compliances remain at
    their healthy values, and Hertz is handled at pair level.
    """
    if crack_depth < 0.0:
        raise ValueError("Crack depth must be non-negative.")
    shear_modulus = context.young / (2.0 * (1.0 + context.poisson))
    h = (gear.root - r_contact) if gear.internal else (r_contact - gear.root)
    h = max(h, 0.05e-3)
    fa, fb = math.sin(alpha_work), math.cos(alpha_work)
    root_thickness = healthy.section_width(gear, gear.root, context.alpha_nom, x_equiv)
    h_c = root_thickness / 2.0
    h_q = h_c - crack_depth * math.sin(crack_angle)

    def half_thickness(x: float) -> float:
        radius = gear.root - x if gear.internal else gear.root + x
        return healthy.section_width(gear, radius, context.alpha_nom, x_equiv) / 2.0

    # The section law is continuous at h_x=h_q but its slope changes there.
    # Providing the physical transition as a quadrature breakpoint prevents
    # round-off warnings at deeper crack cases without smoothing the model.
    transition_points: list[float] = []
    h_at_root, h_at_contact = half_thickness(0.0), half_thickness(h)
    if min(h_at_root, h_at_contact) < h_q < max(h_at_root, h_at_contact):
        transition_points.append(brentq(lambda x: half_thickness(x) - h_q, 0.0, h))

    def sections_at(x: float) -> tuple[float, float, float]:
        radius = gear.root - x if gear.internal else gear.root + x
        healthy_thickness = healthy.section_width(gear, radius, context.alpha_nom, x_equiv)
        h_x = healthy_thickness / 2.0
        if h_x > h_q:
            cracked_thickness = h_c + h_x - crack_depth * math.sin(crack_angle)
            cracked_thickness = max(cracked_thickness, min_ligament_fraction * healthy_thickness)
        else:
            cracked_thickness = healthy_thickness
        area_cracked = face_width * cracked_thickness
        inertia_cracked = face_width * cracked_thickness**3 / 12.0
        area_healthy = face_width * healthy_thickness
        return area_cracked, inertia_cracked, area_healthy

    cb = quad(lambda x: (fb * (h - x) - fa * h) ** 2 / (context.young * sections_at(x)[1]),
              0.0, h, points=transition_points, epsabs=1e-18, epsrel=1e-8, limit=100)[0]
    ca = quad(lambda x: fa**2 / (context.young * sections_at(x)[2]),
              0.0, h, points=transition_points, epsabs=1e-18, epsrel=1e-8, limit=100)[0]
    cs = quad(lambda x: 1.2 * fb**2 / (shear_modulus * sections_at(x)[0]),
              0.0, h, points=transition_points, epsabs=1e-18, epsrel=1e-8, limit=100)[0]

    if foundation_compliance is None:
        lf = max(abs(gear.base - gear.root), 0.05e-3)
        area0 = face_width * root_thickness
        inertia0 = face_width * root_thickness**3 / 12.0
        cf = (fb**2 * lf**3 / (3.0 * context.young * inertia0)
              + fa**2 * lf / (context.young * area0)
              + 1.2 * fb**2 * lf / (shear_modulus * area0))
    else:
        cf = foundation_compliance
    return cb + ca + cs + cf


def pair_stiffness_with_root_crack(
        first: healthy.Gear, second: healthy.Gear, r_first: float, r_second: float,
        face_width: float, damaged_first: bool, damaged_second: bool, crack_depth: float,
        crack_angle: float, alpha_work: float, context: Context, x_first: float = 0.0,
        x_second: float = 0.0, foundation_first: float | None = None,
        foundation_second: float | None = None, min_ligament_fraction: float = 0.05) -> float:
    """Pair stiffness: one cracked tooth is substituted only when it is active."""
    if damaged_first:
        c_first = root_crack_tooth_compliance(first, r_first, face_width, alpha_work, context,
                                               crack_depth, crack_angle, x_first,
                                               foundation_first, min_ligament_fraction)
    else:
        c_first = healthy.tooth_compliance(first, r_first, face_width, alpha_work,
                                           context.alpha_nom, context.young, context.poisson,
                                           x_first, foundation_first)
    if damaged_second:
        c_second = root_crack_tooth_compliance(second, r_second, face_width, alpha_work, context,
                                                crack_depth, crack_angle, x_second,
                                                foundation_second, min_ligament_fraction)
    else:
        c_second = healthy.tooth_compliance(second, r_second, face_width, alpha_work,
                                            context.alpha_nom, context.young, context.poisson,
                                            x_second, foundation_second)
    c_hertz = 1.0 / healthy.hertz_stiffness(context.young, context.poisson, face_width)
    return 1.0 / (c_first + c_second + c_hertz)


def active_pairs(mesh: healthy.MeshPath, display_phase: float, display_cycle: int,
                 phase_origin: float) -> list[tuple[int, float]]:
    """Physical tooth-entry numbers and ages at one displayed sample."""
    raw = display_cycle + display_phase + phase_origin
    entry = math.floor(raw)
    local_phase = raw - entry
    pairs = [(entry, local_phase)]
    if local_phase < mesh.epsilon - 1.0:
        pairs.append((entry - 1, local_phase + 1.0))
    return pairs


def contact_pair(mesh: healthy.MeshPath, first: healthy.Gear, second: healthy.Gear,
                 entry: int, age: float, alpha_work: float) -> tuple[float, float]:
    return healthy.contact_radii(mesh, first, second, age / mesh.epsilon, alpha_work)


def build_local_energy_cache(context: Context, mesh_name: str, crack_depth: float,
                             crack_angle: float, min_ligament_fraction: float,
                             fault_kind: str, planet_ring_offset: int,
                             samples_per_cycle: int) -> list[list[tuple[int, float, float, float]]]:
    """Evaluate every local H/F contact state once, at the requested grid."""
    if mesh_name == "sp":
        mesh, first, second = context.sp_path, context.sun, context.planet
        alpha, face, origin = context.alpha_sp, context.face_sp, context.sp_origin
        first_x, second_x, first_foundation, second_foundation = context.x_each, context.x_each, None, None
    elif mesh_name == "pr":
        mesh, first, second = context.pr_path, context.planet, context.ring
        alpha, face, origin = context.alpha_pr, context.face_pr, context.pr_origin
        first_x, second_x, first_foundation, second_foundation = 0.0, 0.0, None, context.ring_foundation_compliance
    else:
        raise ValueError(mesh_name)
    cache: list[list[tuple[int, float, float, float]]] = []
    for local_sample in range(samples_per_cycle + 1):
        local_display_phase = local_sample / samples_per_cycle
        entries_and_ages = active_pairs(mesh, local_display_phase, 0, origin)
        cached_pairs: list[tuple[int, float, float, float]] = []
        for base_entry, age in entries_and_ages:
            r_first, r_second = contact_pair(mesh, first, second, base_entry, age, alpha)
            if mesh_name == "sp":
                damage_first, damage_second = fault_kind == "sun", fault_kind == "planet"
            else:
                damage_first, damage_second = fault_kind == "planet", False
            # Use the existing validated healthy evaluator directly rather
            # than recomputing an equivalent zero-depth crack case.
            healthy_pair = healthy.pair_stiffness(
                first, second, r_first, r_second, face, alpha, context.alpha_nom,
                context.young, context.poisson, first_x, second_x,
                first_foundation, second_foundation)
            damaged_pair = pair_stiffness_with_root_crack(
                first, second, r_first, r_second, face, damage_first, damage_second, crack_depth, crack_angle,
                alpha, context, first_x, second_x, first_foundation, second_foundation, min_ligament_fraction)
            cached_pairs.append((base_entry, age, healthy_pair, damaged_pair))
        cache.append(cached_pairs)
    return cache


def tooth_number(entry: int, teeth: int, offset: int = 0, sign: int = 1) -> int:
    """One-based tooth identity; sign=-1 follows the planet rotation sense."""
    return (sign * entry + offset) % teeth + 1


def make_series(context: Context, mesh_name: str, cycles: int, crack_depth: float,
                crack_angle: float, min_ligament_fraction: float, fault_kind: str,
                planet_index: int, planet_ring_offset: int,
                samples_per_cycle: int,
                local_cache: list[list[tuple[int, float, float, float]]] | None = None) -> list[dict[str, float | int]]:
    """Create a local-fault and matched-healthy series over integer recurrence cycles."""
    if mesh_name == "sp":
        mesh, first, second = context.sp_path, context.sun, context.planet
        alpha, face, origin = context.alpha_sp, context.face_sp, context.sp_origin
        first_x, second_x = context.x_each, context.x_each
        first_foundation, second_foundation = None, None
        planet_offset = 0
        sun_offset = (planet_index - 1) * (context.zs // context.n_planet)
    elif mesh_name == "pr":
        mesh, first, second = context.pr_path, context.planet, context.ring
        alpha, face, origin = context.alpha_pr, context.face_pr, context.pr_origin
        first_x, second_x = 0.0, 0.0
        first_foundation, second_foundation = None, context.ring_foundation_compliance
        planet_offset = planet_ring_offset
        sun_offset = 0
    else:
        raise ValueError(mesh_name)

    # A contact-radius/energy result depends on local phase and on whether
    # the tooth is damaged, not on its absolute tooth number.  Cache the
    # exactly evaluated 2001 local samples, then select H/F by the real tooth
    # index in each recurrence cycle.  This is neither interpolation nor a
    # repeated/fitted waveform.
    if local_cache is None:
        local_cache = build_local_energy_cache(context, mesh_name, crack_depth, crack_angle,
                                                min_ligament_fraction, fault_kind,
                                                planet_ring_offset, samples_per_cycle)

    rows: list[dict[str, float | int]] = []
    for sample in range(cycles * samples_per_cycle + 1):
        display_total = sample / samples_per_cycle
        display_cycle = min(int(math.floor(display_total)), cycles - 1)
        display_phase = display_total - math.floor(display_total)
        # At the final endpoint retain the physical display phase=1 within the
        # preceding displayed cycle.  This gives an unambiguous endpoint and
        # preserves exact healthy-cycle closure for a complete recurrence.
        if sample == cycles * samples_per_cycle:
            display_cycle, display_phase = cycles - 1, 1.0
        cached_pairs = local_cache[int(round(display_phase * samples_per_cycle))]
        pair_h: list[float] = []
        pair_f: list[float] = []
        fault_flags: list[int] = []
        entries: list[int] = []
        first_teeth: list[int] = []
        second_teeth: list[int] = []
        for base_entry, age, healthy_pair, damaged_pair in cached_pairs:
            entry = display_cycle + base_entry
            if mesh_name == "sp":
                first_tooth = tooth_number(entry, context.zs, sun_offset, sign=1)
                second_tooth = tooth_number(entry, context.zp, 0, sign=-1)
                damaged_first = fault_kind == "sun" and first_tooth == 1
                damaged_second = fault_kind == "planet" and second_tooth == 1
            else:
                first_tooth = tooth_number(entry, context.zp, planet_offset, sign=-1)
                second_tooth = tooth_number(entry, context.zr, (planet_index - 1) * (context.zr // context.n_planet), sign=1)
                damaged_first = fault_kind == "planet" and first_tooth == 1
                damaged_second = False
            faulty_pair = damaged_pair if damaged_first or damaged_second else healthy_pair
            pair_h.append(healthy_pair)
            pair_f.append(faulty_pair)
            fault_flags.append(int(damaged_first or damaged_second))
            entries.append(entry)
            first_teeth.append(first_tooth)
            second_teeth.append(second_tooth)
        row: dict[str, float | int] = {
            "mesh_phase_cycles": display_total,
            "mesh_cycle_index": int(math.floor(display_total)) if sample < cycles * samples_per_cycle else cycles,
            "phase_norm": display_phase,
            "active_pair_count": len(cached_pairs),
            "fault_active_pair_count": sum(fault_flags),
            "k_healthy_N_per_m": sum(pair_h),
            "k_fault_N_per_m": sum(pair_f),
            "delta_k_N_per_m": sum(pair_h) - sum(pair_f),
        }
        for number in (0, 1):
            suffix = number + 1
            if number < len(cached_pairs):
                row.update({
                    f"entry_pair_{suffix}": entries[number],
                    f"age_pair_{suffix}": cached_pairs[number][1],
                    f"first_tooth_pair_{suffix}": first_teeth[number],
                    f"second_tooth_pair_{suffix}": second_teeth[number],
                    f"is_fault_pair_{suffix}": fault_flags[number],
                    f"k_pair_healthy_{suffix}_N_per_m": pair_h[number],
                    f"k_pair_fault_{suffix}_N_per_m": pair_f[number],
                })
            else:
                row.update({
                    f"entry_pair_{suffix}": "", f"age_pair_{suffix}": "",
                    f"first_tooth_pair_{suffix}": "", f"second_tooth_pair_{suffix}": "",
                    f"is_fault_pair_{suffix}": "", f"k_pair_healthy_{suffix}_N_per_m": "",
                    f"k_pair_fault_{suffix}_N_per_m": "",
                })
        rows.append(row)
    return rows


def write_rows(path: Path, rows: list[dict[str, float | int]]) -> None:
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def numeric(rows: list[dict[str, float | int]], field: str) -> np.ndarray:
    return np.asarray([float(row[field]) for row in rows])


def configure_plot() -> None:
    mpl.rcParams.update({
        "font.family": "Times New Roman", "font.size": 9.5,
        "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "legend.frameon": False,
    })


def apply_tick_font(ax: plt.Axes) -> None:
    for text in (*ax.get_xticklabels(), *ax.get_yticklabels()):
        text.set_fontproperties(TIMES)


def set_mixed_ylabel(ax: plt.Axes, chinese_label: str) -> None:
    """Render Chinese separately; draw all unit glyphs directly with times.ttf."""
    ax.set_ylabel("")
    common = {"transform": ax.transAxes, "rotation": 90, "rotation_mode": "anchor",
              "ha": "center", "va": "center", "clip_on": False}
    ax.text(-0.075, 0.34, chinese_label, fontproperties=SIMSUN, **common)
    ax.text(-0.075, 0.73, "(10⁸ N/m)", fontproperties=TIMES, **common)


def plot_sun(out: Path, q: float, by_planet: dict[int, list[dict[str, float | int]]]) -> None:
    """Quantitative-grid figure: three SP meshes demonstrate localized loss."""
    fig, axes = plt.subplots(2, 1, figsize=(6.1, 4.8), dpi=600, sharex=True,
                             gridspec_kw={"height_ratios": [2.0, 1.0]})
    colors = ["#1f4e79", "#4f81bd", "#79a9d1"]
    for planet, color in zip(sorted(by_planet), colors):
        rows = by_planet[planet]
        x = numeric(rows, "mesh_phase_cycles")
        healthy_k = numeric(rows, "k_healthy_N_per_m") / 1e8
        faulty_k = numeric(rows, "k_fault_N_per_m") / 1e8
        delta = numeric(rows, "delta_k_N_per_m") / 1e8
        if planet == 1:
            axes[0].plot(x, healthy_k, color="#4b4b4b", linewidth=0.8, label="健康")
            axes[0].plot(x, faulty_k, color=color, linewidth=1.05, label="故障 P1 网格")
        else:
            axes[0].plot(x, faulty_k, color=color, linewidth=0.85, label=f"故障 P{planet} 网格")
        axes[1].plot(x, delta, color=color, linewidth=0.95, label=f"P{planet}")
    set_mixed_ylabel(axes[0], "SP 啮合刚度")
    set_mixed_ylabel(axes[1], "刚度损失")
    axes[1].set_xlabel("显示啮合相位（周期）", fontproperties=SIMSUN)
    axes[0].legend(prop=SIMSUN, loc="upper right", ncol=2, fontsize=8)
    axes[1].legend(prop=TIMES, loc="upper right", ncol=3, fontsize=8)
    for ax in axes:
        ax.margins(x=0)
        apply_tick_font(ax)
    fig.tight_layout(pad=0.65, h_pad=0.35)
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)


def plot_planet(out: Path, q: float, sp: list[dict[str, float | int]],
                pr: list[dict[str, float | int]]) -> None:
    """Two-panel global view retaining the SP/PR fault-tooth phase relation."""
    fig, axes = plt.subplots(2, 1, figsize=(6.1, 4.8), dpi=600, sharex=True)
    for ax, rows, label, color in zip(axes, (sp, pr), ("SP", "PR"), ("#1f4e79", "#a64b00")):
        x = numeric(rows, "mesh_phase_cycles")
        h = numeric(rows, "k_healthy_N_per_m") / 1e8
        f = numeric(rows, "k_fault_N_per_m") / 1e8
        ax.plot(x, h, color="#4b4b4b", linewidth=0.75, label="健康")
        ax.plot(x, f, color=color, linewidth=1.0, label="故障")
        set_mixed_ylabel(ax, f"{label} 啮合刚度")
        ax.legend(prop=SIMSUN, loc="upper right", fontsize=8)
        ax.margins(x=0)
        apply_tick_font(ax)
    axes[1].set_xlabel("显示啮合相位（周期）", fontproperties=SIMSUN)
    fig.tight_layout(pad=0.65, h_pad=0.45)
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)


def reference_csv_check(context: Context) -> float:
    """Ensure the healthy SP implementation remains byte-for-value compatible."""
    phase, _, _, expected, _ = healthy.mesh_curve(
        context.sp_path, context.sun, context.planet, context.face_sp, context.alpha_sp,
        context.alpha_nom, context.young, context.poisson, context.x_each, context.x_each,
        context.sp_origin,
    )
    with (TVMS_OUT / "tvms_sp.csv").open(newline="", encoding="utf-8") as handle:
        saved = np.asarray([float(row["k_total_N_per_m"]) for row in csv.DictReader(handle)])
    if not np.allclose(expected, saved, rtol=2e-11, atol=0.0):
        raise RuntimeError("Healthy SP implementation no longer matches the saved validated TVMS CSV.")
    return float(np.max(np.abs(expected - saved)))


def check_series(rows: list[dict[str, float | int]], crack_depth: float, label: str) -> list[str]:
    h, f, delta = (numeric(rows, key) for key in ("k_healthy_N_per_m", "k_fault_N_per_m", "delta_k_N_per_m"))
    flags = np.asarray([int(row["fault_active_pair_count"]) for row in rows])
    if not (np.all(np.isfinite(h)) and np.all(np.isfinite(f)) and np.all(np.isfinite(delta))):
        raise RuntimeError(f"{label}: NaN/Inf found")
    if np.any(h <= 0) or np.any(f < 0) or np.any(delta < -1e-5):
        raise RuntimeError(f"{label}: nonphysical stiffness/loss found")
    outside = delta[flags == 0]
    if outside.size and not np.allclose(outside, 0.0, atol=1e-7, rtol=0.0):
        raise RuntimeError(f"{label}: stiffness changed outside a fault-tooth contact")
    if crack_depth > 0 and (not np.any(flags > 0) or not np.any(delta[flags > 0] > 0)):
        raise RuntimeError(f"{label}: fault tooth generated no local stiffness loss")
    return [
        f"{label}: finite values=PASS; local-only loss=PASS; min healthy/fault={h.min():.6e}/{f.min():.6e} N/m; max delta={delta.max():.6e} N/m.",
    ]


def fault_entries(rows: list[dict[str, float | int]]) -> list[int]:
    result: set[int] = set()
    for row in rows:
        for number in (1, 2):
            if row[f"is_fault_pair_{number}"] == 1:
                result.add(int(row[f"entry_pair_{number}"]))
    return sorted(result)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth-mm", type=float, action="append", default=None,
                        help="Root-crack extension q_c in mm; may be repeated.")
    parser.add_argument("--run-name", default=None, help="Unique results subdirectory name.")
    args = parser.parse_args()
    config = read_config()
    depth_values_mm = config["crack_depth_mm_values"] if args.depth_mm is None else args.depth_mm
    if any(depth < 0.0 for depth in depth_values_mm):
        raise ValueError("Every crack depth must be non-negative.")
    if config["damage_model"] != "straight_root_crack_energy":
        raise RuntimeError("This runner only implements the configured straight root-crack energy model.")
    crack_angle = math.radians(float(config["crack_angle_deg"]))
    min_ligament_fraction = float(config["min_ligament_fraction"])
    run_name = args.run_name or f"crack_tvms_{datetime.now():%Y%m%d_%H%M%S}"
    out = TVMS_OUT / run_name
    if out.exists():
        raise RuntimeError(f"Refusing to overwrite existing result directory: {out}")
    out.mkdir(parents=True)
    context = build_context()
    healthy_difference = reference_csv_check(context)
    samples = int(config["samples_per_mesh_cycle"])
    fault_tooth, target_planet, pr_offset = int(config["fault_tooth"]), int(config["target_planet"]), int(config["planet_ring_tooth_offset"])
    if fault_tooth != 1:
        raise RuntimeError("Current tooth-number mapping is explicitly validated for fault tooth 1 only.")
    configure_plot()
    report = [
        "Parameterized local straight-root-crack potential-energy TVMS", "",
        "Scope: TVMS calculation only; no dynamic response calculation was run.",
        "Fault form: straight root crack on one tooth, following Peng Yue (2025), Eqs. (3-1)--(3-3).",
        "q_c: crack extension from the tooth root; gamma: angle between crack and tooth centreline.",
        "Recomputed on the named fault tooth only: bending and shear compliances from crack-reduced I_x and A_x.",
        "Held at validated healthy values: axial compression, tooth foundation, mating tooth, Hertz contact term, contact geometry, active-pair logic, contact ratio, working pressure angle and ring FEM foundation compliance.",
        f"gamma={math.degrees(crack_angle):.6f} deg; minimum retained ligament fraction={min_ligament_fraction:.6f}.",
        f"SP alpha={math.degrees(context.alpha_sp):.9f} deg; PR alpha={math.degrees(context.alpha_pr):.9f} deg; epsilon_SP={context.sp_path.epsilon:.9f}; epsilon_PR={context.pr_path.epsilon:.9f}.",
        f"B_SP={context.face_sp:.9e} m; B_PR={context.face_pr:.9e} m; C_f,ring={context.ring_foundation_compliance:.9e} m/N.",
        f"Saved healthy SP reference agreement: max absolute difference={healthy_difference:.6e} N/m (PASS).",
        f"Sun fault recurrence={context.zs} displayed mesh cycles; planet fault recurrence={context.zp} displayed mesh cycles; PR planet-tooth entry offset={pr_offset}.",
        "",
    ]
    for depth_mm in depth_values_mm:
        crack_depth = float(depth_mm) * 1e-3
        tag = f"qc{str(f'{depth_mm:.3f}').replace('.', 'p')}mm"
        sun_sp_cache = build_local_energy_cache(context, "sp", crack_depth, crack_angle, min_ligament_fraction,
                                                 "sun", pr_offset, samples)
        planet_sp_cache = build_local_energy_cache(context, "sp", crack_depth, crack_angle, min_ligament_fraction,
                                                    "planet", pr_offset, samples)
        planet_pr_cache = build_local_energy_cache(context, "pr", crack_depth, crack_angle, min_ligament_fraction,
                                                    "planet", pr_offset, samples)
        sun_by_planet: dict[int, list[dict[str, float | int]]] = {}
        for planet_index in range(1, context.n_planet + 1):
            series = make_series(context, "sp", context.zs, crack_depth, crack_angle, min_ligament_fraction,
                                 "sun", planet_index, pr_offset, samples, sun_sp_cache)
            sun_by_planet[planet_index] = series
            write_rows(out / f"sun_fault_sp_p{planet_index}_{tag}.csv", series)
            report.extend(check_series(series, crack_depth, f"q_c={depth_mm:g} mm, sun-tooth crack, SP P{planet_index}"))
        planet_sp = make_series(context, "sp", context.zp, crack_depth, crack_angle, min_ligament_fraction,
                                "planet", target_planet, pr_offset, samples, planet_sp_cache)
        planet_pr = make_series(context, "pr", context.zp, crack_depth, crack_angle, min_ligament_fraction,
                                "planet", target_planet, pr_offset, samples, planet_pr_cache)
        write_rows(out / f"planet_fault_sp_p{target_planet}_{tag}.csv", planet_sp)
        write_rows(out / f"planet_fault_pr_p{target_planet}_{tag}.csv", planet_pr)
        report.extend(check_series(planet_sp, crack_depth, f"q_c={depth_mm:g} mm, planet-tooth crack, SP P{target_planet}"))
        report.extend(check_series(planet_pr, crack_depth, f"q_c={depth_mm:g} mm, planet-tooth crack, PR P{target_planet}"))
        sp_entries, pr_entries = fault_entries(planet_sp), fault_entries(planet_pr)
        if not sp_entries or not pr_entries or any((pr_entry - sp_entries[0]) % context.zp != pr_offset % context.zp for pr_entry in pr_entries):
            raise RuntimeError("Planet SP/PR tooth-phase relation check failed.")
        report += [
            f"q_c={depth_mm:g} mm: planet tooth-1 physical entry indices: SP={sp_entries}; PR={pr_entries}; offset modulo Zp={(pr_entries[0]-sp_entries[0]) % context.zp} (PASS).",
            f"q_c={depth_mm:g} mm: depth/module ratio={crack_depth / (context.sun.pitch * 2 / context.zs):.6f}; gamma={math.degrees(crack_angle):.6f} deg.",
        ]
        plot_sun(out / f"Fig_crack_sun_{tag}.png", crack_depth, sun_by_planet)
        plot_planet(out / f"Fig_crack_planet_{tag}.png", crack_depth, planet_sp, planet_pr)
    (out / "fault_tvms_summary.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))
    print(f"Output directory: {out}")


if __name__ == "__main__":
    main()
