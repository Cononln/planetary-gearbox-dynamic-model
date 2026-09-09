"""Healthy spur-gear TVMS using a geometry-resolved potential-energy model.

This runner is deliberately separate from the existing mean-stiffness model.
It reads frozen values in ``src/pg_parameters.m`` and does not alter the
dynamic model. The user-authorised 39.75 mm SP working centre is represented
as an equivalent tooth-thickness correction, not a measured profile shift.
"""
from __future__ import annotations

import csv
import ast
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "tvms"
OUT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".matplotlib"))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from scipy.integrate import quad
from ring_foundation_fem import RingRimResult, annular_rim_normal_compliance

PFILE = ROOT / "src" / "pg_parameters.m"
N_PHASE = 2001


def source_value(name: str) -> float:
    """Read a numeric SI-valued assignment p.<name> from MATLAB."""
    source = PFILE.read_text(encoding="utf-8")
    match = re.search(
        rf"^\s*p\.{re.escape(name)}\s*=\s*(.+?)\s*;",
        source,
        flags=re.MULTILINE,
    )
    if not match:
        raise RuntimeError(f"Missing literal parameter p.{name} in {PFILE}")
    expression = match.group(1).strip()
    expression = re.sub(r"deg2rad\(([-+]?\d+(?:\.\d*)?)\)",
                        lambda item: str(math.radians(float(item.group(1)))), expression)
    tree = ast.parse(expression, mode="eval")
    allowed = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Add, ast.Sub,
               ast.Mult, ast.Div, ast.USub, ast.UAdd, ast.Constant)
    if not all(isinstance(node, allowed) for node in ast.walk(tree)):
        raise RuntimeError(f"Unsupported non-numeric parameter expression for p.{name}: {expression}")
    return float(eval(compile(tree, "<parameter>", "eval"), {"__builtins__": {}}, {}))


def involute(alpha: float) -> float:
    return math.tan(alpha) - alpha


@dataclass(frozen=True)
class Gear:
    name: str
    z: int
    pitch: float
    base: float
    tip: float
    root: float
    internal: bool = False


def section_width(gear: Gear, r: float, alpha_nom: float, x_equiv: float = 0.0) -> float:
    """Involute tooth thickness; internal gear uses the opposite sign."""
    r_eff = max(r, gear.base * (1.0 + 1e-12))
    alpha_r = math.acos(min(1.0, gear.base / r_eff))
    sign = -1.0 if gear.internal else 1.0
    half_angle = math.pi / (2.0 * gear.z) + sign * (involute(alpha_nom) - involute(alpha_r))
    half_angle += x_equiv * math.tan(alpha_nom) / gear.z
    return max(0.05e-3, 2.0 * r_eff * abs(half_angle))


def tooth_compliance(gear: Gear, r_contact: float, face_width: float,
                     alpha_work: float, alpha_nom: float, young: float,
                     poisson: float, x_equiv: float = 0.0,
                     foundation_compliance: float | None = None) -> float:
    """Unit-load bending, axial, shear and root-section compliances.

    Bending uses [Fb(d-x)-Fa*h]^2/(E I_x). The foundation term is a
    disclosed root-section energy continuation because no fillet/rim model
    is available.
    """
    shear_modulus = young / (2.0 * (1.0 + poisson))
    h = (gear.root - r_contact) if gear.internal else (r_contact - gear.root)
    h = max(h, 0.05e-3)
    fa, fb = math.sin(alpha_work), math.cos(alpha_work)

    def section_at(x: float) -> tuple[float, float]:
        radius = gear.root - x if gear.internal else gear.root + x
        thickness = section_width(gear, radius, alpha_nom, x_equiv)
        return face_width * thickness, face_width * thickness**3 / 12.0

    cb = quad(lambda x: (fb * (h - x) - fa * h) ** 2 / (young * section_at(x)[1]),
              0.0, h, epsabs=1e-18, epsrel=1e-8, limit=100)[0]
    ca = quad(lambda x: fa**2 / (young * section_at(x)[0]),
              0.0, h, epsabs=1e-18, epsrel=1e-8, limit=100)[0]
    cs = quad(lambda x: 1.2 * fb**2 / (shear_modulus * section_at(x)[0]),
              0.0, h, epsabs=1e-18, epsrel=1e-8, limit=100)[0]
    if foundation_compliance is None:
        area0, inertia0 = section_at(0.0)
        lf = max(abs(gear.base - gear.root), 0.05e-3)
        cf = (fb**2 * lf**3 / (3.0 * young * inertia0)
              + fa**2 * lf / (young * area0)
              + 1.2 * fb**2 * lf / (shear_modulus * area0))
    else:
        cf = foundation_compliance
    return cb + ca + cs + cf


def hertz_stiffness(young: float, poisson: float, face_width: float) -> float:
    """Line-contact Hertz stiffness (N/m) under the unit-load convention."""
    equivalent_modulus = 1.0 / (2.0 * (1.0 - poisson**2) / young)
    return math.pi * face_width * equivalent_modulus


def pair_stiffness(first: Gear, second: Gear, r_first: float, r_second: float,
                   face_width: float, alpha_work: float, alpha_nom: float,
                   young: float, poisson: float, x_first: float = 0.0,
                   x_second: float = 0.0, foundation_first: float | None = None,
                   foundation_second: float | None = None) -> float:
    c1 = tooth_compliance(first, r_first, face_width, alpha_work, alpha_nom, young, poisson, x_first, foundation_first)
    c2 = tooth_compliance(second, r_second, face_width, alpha_work, alpha_nom, young, poisson, x_second, foundation_second)
    return 1.0 / (c1 + c2 + 1.0 / hertz_stiffness(young, poisson, face_width))


@dataclass(frozen=True)
class MeshPath:
    kind: str
    epsilon: float
    length: float
    start: float
    end: float


def external_path(sun: Gear, planet: Gear, alpha_work: float, base_pitch: float) -> MeshPath:
    """External-mesh path of contact, referenced to the pitch point."""
    sun_recess = math.sqrt(sun.tip**2 - sun.base**2) - sun.base * math.tan(alpha_work)
    planet_approach = math.sqrt(planet.tip**2 - planet.base**2) - planet.base * math.tan(alpha_work)
    length = sun_recess + planet_approach
    return MeshPath("external", length / base_pitch, length, -planet_approach, sun_recess)


def internal_path(planet: Gear, ring: Gear, alpha_work: float, centre: float,
                  base_pitch: float) -> MeshPath:
    """Internal-mesh contact length; this is not the external formula."""
    length = (centre * math.sin(alpha_work)
              - math.sqrt(ring.tip**2 - ring.base**2)
              + math.sqrt(planet.tip**2 - planet.base**2))
    return MeshPath("internal", length / base_pitch, length, 0.0, length)


def contact_radii(mesh: MeshPath, first: Gear, second: Gear, fraction: float,
                  alpha_work: float) -> tuple[float, float]:
    """External and internal gears receive separate contact trajectories."""
    fraction = min(1.0, max(0.0, fraction))
    if mesh.kind == "external":
        s = mesh.start + fraction * mesh.length
        l1 = first.base * math.tan(alpha_work) + s
        l2 = second.base * math.tan(alpha_work) - s
        r1 = math.sqrt(first.base**2 + max(0.0, l1)**2)
        r2 = math.sqrt(second.base**2 + max(0.0, l2)**2)
    else:
        # Same-side internal-gear involute motion: planet moves towards its
        # tip, whereas the ring contact moves from its tip towards its root.
        r1 = first.base + fraction * (first.tip - first.base)
        r2 = second.tip + fraction * (second.root - second.tip)
    r1 = min(max(r1, first.base), first.tip)
    r2 = min(max(r2, second.tip if second.internal else second.base),
             second.root if second.internal else second.tip)
    return r1, r2


def mesh_curve(mesh: MeshPath, first: Gear, second: Gear, face_width: float,
               alpha_work: float, alpha_nom: float, young: float, poisson: float,
               x_first: float = 0.0, x_second: float = 0.0,
               phase_origin: float = 0.0, foundation_first: float | None = None,
               foundation_second: float | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate one mesh period after an exact cyclic phase-origin change.

    ``phase_origin`` only changes which physical mesh phase is reported as
    zero.  No stiffness interpolation, smoothing, or load-transfer model is
    introduced here.
    """
    phase = np.linspace(0.0, 1.0, N_PHASE)
    stiffness = np.empty_like(phase)
    pair_1 = np.full_like(phase, np.nan)
    pair_2 = np.full_like(phase, np.nan)
    active_count = np.zeros_like(phase, dtype=int)
    for idx, display_phi in enumerate(phase):
        phi = (display_phi + phase_origin) % 1.0
        ages = [phi]  # tooth pair entering in the current mesh period
        if phi < mesh.epsilon - 1.0:
            ages.append(phi + 1.0)  # preceding tooth pair still in contact
        pairs = []
        for age in ages:
            r1, r2 = contact_radii(mesh, first, second, age / mesh.epsilon, alpha_work)
            pairs.append(pair_stiffness(first, second, r1, r2, face_width, alpha_work,
                                        alpha_nom, young, poisson, x_first, x_second,
                                        foundation_first, foundation_second))
        pair_1[idx] = pairs[0]
        if len(pairs) == 2:
            pair_2[idx] = pairs[1]
        active_count[idx] = len(pairs)
        stiffness[idx] = sum(pairs)
    return phase, pair_1, pair_2, stiffness, active_count


def configure_matplotlib() -> None:
    try:
        font_manager.findfont("Times New Roman", fallback_to_default=False)
    except ValueError as exc:
        raise RuntimeError("Times New Roman is required for the requested thesis figures but is not installed.") from exc
    mpl.rcParams.update({"font.family": "Times New Roman", "font.size": 10,
                         "axes.linewidth": 0.8, "axes.spines.top": False,
                         "axes.spines.right": False, "xtick.direction": "out",
                         "ytick.direction": "out", "legend.frameon": False,
                         "savefig.facecolor": "white", "figure.facecolor": "white"})


def shade_double_pair_regions(ax: plt.Axes, phase: np.ndarray, active_count: np.ndarray) -> None:
    """Mark computed two-pair intervals; this only annotates the result."""
    in_region = False
    start = 0.0
    for index, is_double in enumerate(active_count == 2):
        if is_double and not in_region:
            start = phase[index]
            in_region = True
        if in_region and (not is_double or index == len(phase) - 1):
            stop = phase[index] if not is_double else phase[index]
            ax.axvspan(start, stop, facecolor="#e8eef5", edgecolor="none", zorder=0)
            ax.axvline(start, color="#9aa7b5", linestyle="--", linewidth=0.65, zorder=1)
            ax.axvline(stop, color="#9aa7b5", linestyle="--", linewidth=0.65, zorder=1)
            in_region = False


def save_cycle_figure(filename: str, phase: np.ndarray, pair_1: np.ndarray,
                      pair_2: np.ndarray, total: np.ndarray,
                      active_count: np.ndarray, mesh_label: str) -> None:
    """Formal one-period thesis plot: pair terms and their direct sum."""
    fig, ax = plt.subplots(figsize=(5.4, 3.45), dpi=330)
    shade_double_pair_regions(ax, phase, active_count)
    ax.plot(phase, pair_1 / 1e8, linewidth=0.95, color="#7c8794", label="Active pair 1")
    ax.plot(phase, pair_2 / 1e8, linewidth=0.95, color="#aab4bf", label="Active pair 2")
    ax.plot(phase, total / 1e8, linewidth=1.45, color="#1f4e79", label=f"Total {mesh_label}")
    ax.set(xlim=(0.0, 1.0), xticks=np.arange(0.0, 1.1, 0.2),
           xlabel="Normalized meshing phase",
           ylabel=r"Meshing stiffness ($10^8$ N/m)")
    ax.legend(loc="best", fontsize=8.2)
    fig.tight_layout(pad=0.55)
    fig.savefig(OUT / filename, dpi=330, bbox_inches="tight")
    plt.close(fig)


def save_two_cycle_figure(filename: str, phase: np.ndarray, sp: np.ndarray,
                          pr: np.ndarray) -> None:
    """Two exact adjacent periods, retained as a periodicity diagnostic."""
    fig, axes = plt.subplots(2, 1, figsize=(5.4, 4.7), dpi=330, sharex=True)
    for ax, signal, label, color in zip(axes, (sp, pr), ("Sun--planet", "Planet--ring"), ("#1f4e79", "#a64b00")):
        ax.plot(phase, signal / 1e8, color=color, linewidth=1.2)
        ax.plot(phase + 1.0, signal / 1e8, color=color, linewidth=1.2)
        ax.axvline(1.0, color="#9aa7b5", linestyle="--", linewidth=0.65)
        ax.set_ylabel(rf"{label}\n($10^8$ N/m)")
    axes[-1].set(xlim=(0.0, 2.0), xticks=np.arange(0.0, 2.1, 0.5),
                 xlabel="Normalized meshing phase (two periods)")
    fig.tight_layout(pad=0.6)
    fig.savefig(OUT / filename, dpi=330, bbox_inches="tight")
    plt.close(fig)


def save_extended_total_figure(filename: str, phase: np.ndarray, total: np.ndarray,
                               mesh_label: str, color: str, end_phase: float = 1.65) -> None:
    """Formal TVMS figure from exact periodic extension, k(psi+1)=k(psi).

    The extension only concatenates already computed samples. It uses no
    interpolation, filtering, load-transfer weighting, or recomputation of
    tooth-pair stiffness.
    """
    x_blocks, y_blocks = [], []
    for period in range(math.ceil(end_phase)):
        mask = phase + period <= end_phase + 1e-12
        x_blocks.append(phase[mask] + period)
        y_blocks.append(total[mask])
    x_plot, y_plot = np.concatenate(x_blocks), np.concatenate(y_blocks)
    fig, ax = plt.subplots(figsize=(5.4, 3.45), dpi=330)
    ax.plot(x_plot, y_plot / 1e8, linewidth=1.35, color=color)
    ax.set(xlim=(0.0, end_phase), xticks=[0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.65],
           xlabel="Normalized meshing phase", ylabel=r"Meshing stiffness ($10^8$ N/m)")
    fig.tight_layout(pad=0.55)
    fig.savefig(OUT / filename, dpi=330, bbox_inches="tight")
    plt.close(fig)


def save_total_figure(filename: str, phase: np.ndarray,
                      curves: list[tuple[np.ndarray, str, str]]) -> None:
    """One-period comparison figure for total mesh stiffness only."""
    fig, ax = plt.subplots(figsize=(5.4, 3.45), dpi=330)
    for signal, label, color in curves:
        ax.plot(phase, signal / 1e8, linewidth=1.25, color=color, label=label)
    ax.set(xlim=(0.0, 1.0), xticks=np.arange(0.0, 1.1, 0.2),
           xlabel="Normalized meshing phase", ylabel=r"Meshing stiffness ($10^8$ N/m)")
    ax.legend(loc="best", fontsize=8.2)
    fig.tight_layout(pad=0.55)
    fig.savefig(OUT / filename, dpi=330, bbox_inches="tight")
    plt.close(fig)


def count_sequence(active_count: np.ndarray) -> list[int]:
    return [int(active_count[0]), *[int(active_count[index]) for index in range(1, len(active_count)) if active_count[index] != active_count[index - 1]]]


def double_interval_text(phase: np.ndarray, active_count: np.ndarray) -> str:
    points = phase[active_count == 2]
    if len(points) == 0:
        return "none"
    return f"[{points[0]:.9f}, {points[-1]:.9f}]"


def write_csv(filename: str, phase: np.ndarray, pair_1: np.ndarray,
              pair_2: np.ndarray, total: np.ndarray, active_count: np.ndarray,
              mesh_label: str, fmesh: float) -> None:
    with (OUT / filename).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["phase_norm", "phase_rad", "time_s", "active_pair_count",
                         "k_pair_1_N_per_m", "k_pair_2_N_per_m", "k_total_N_per_m",
                         f"k_{mesh_label}_N_per_m"])
        for u, first, second, k, count in zip(phase, pair_1, pair_2, total, active_count):
            writer.writerow([f"{u:.12g}", f"{2 * math.pi * u:.12g}",
                             f"{u / fmesh:.12g}", int(count), f"{first:.12g}",
                             "" if np.isnan(second) else f"{second:.12g}",
                             f"{k:.12g}", f"{k:.12g}"])


def stats(signal: np.ndarray) -> str:
    return f"min={signal.min():.6e} N/m, max={signal.max():.6e} N/m, mean={signal.mean():.6e} N/m"


def main() -> None:
    zs, zp, zr = (int(source_value(item)) for item in ("gear.zs", "gear.zp", "gear.zr"))
    module, alpha_nom, centre = (source_value(item) for item in ("gear.module", "gear.alpha_nominal", "gear.carrierRadius"))
    face_sp, face_pr = (source_value(item) for item in ("gear.faceWidthSunPlanet", "gear.faceWidthRing"))
    young, poisson, density = (source_value(item) for item in ("material.gear.E", "material.gear.nu", "material.gear.rho"))
    sun_rpm, n_planet = source_value("operating.sunRpm"), int(source_value("model.nPlanet"))
    sun = Gear("sun", zs, module * zs / 2, module * zs * math.cos(alpha_nom) / 2,
               source_value("gear.addendumRadiusSun"), source_value("gear.rootRadiusSun"))
    planet = Gear("planet", zp, module * zp / 2, module * zp * math.cos(alpha_nom) / 2,
                  source_value("gear.addendumRadiusPlanet"), source_value("gear.rootRadiusPlanet"))
    ring = Gear("ring", zr, module * zr / 2, module * zr * math.cos(alpha_nom) / 2,
                source_value("gear.addendumRadiusRing"), source_value("gear.rootRadiusRing"), True)
    alpha_sp = math.acos((sun.pitch + planet.pitch) * math.cos(alpha_nom) / centre)
    alpha_pr = math.acos((ring.pitch - planet.pitch) * math.cos(alpha_nom) / centre)
    base_pitch = math.pi * module * math.cos(alpha_nom)
    sp_path, pr_path = external_path(sun, planet, alpha_sp, base_pitch), internal_path(planet, ring, alpha_pr, centre, base_pitch)
    if min(sp_path.epsilon, pr_path.epsilon) <= 1.0:
        raise RuntimeError(f"Contact ratio must exceed one: SP={sp_path.epsilon:.9f}, PR={pr_path.epsilon:.9f}")
    x_total = (zs + zp) * (involute(alpha_sp) - involute(alpha_nom)) / (2 * math.tan(alpha_nom))
    x_each = x_total / 2
    # New user-supplied ring geometry supports an annular-rim plane-stress
    # foundation calculation.  The 8 through holes are removed explicitly;
    # the housing-integral outer rim is fixed.  Fine/coarse meshes provide a
    # numerical convergence check without altering the other TVMS terms.
    ring_root_thickness = section_width(ring, ring.root, alpha_nom)
    rim_args = dict(root_radius=ring.root, outer_radius=source_value("gear.outerRadiusRing"),
                    face_width=face_pr, young=young, poisson=poisson,
                    hole_radius=source_value("gear.ringHoleRadius"),
                    hole_diameter=source_value("gear.ringHoleDiameter"),
                    hole_count=int(source_value("gear.ringHoleCount")), tooth_count=zr,
                    root_tooth_thickness=ring_root_thickness, pressure_angle=alpha_pr)
    rim_coarse: RingRimResult = annular_rim_normal_compliance(**rim_args, nr=28, ntheta=504)
    rim_fine: RingRimResult = annular_rim_normal_compliance(**rim_args, nr=36, ntheta=672)
    rim_convergence = abs(rim_fine.compliance_mean-rim_coarse.compliance_mean)/rim_fine.compliance_mean
    if rim_convergence > .05:
        raise RuntimeError(f"Ring-rim FEM mesh convergence failed: {rim_convergence:.3%}")
    ring_foundation_compliance = rim_fine.compliance_mean
    # For 1 < epsilon < 2, origin=epsilon/2 lies in the single-pair region.
    # The displayed interval therefore contains single -> double -> single,
    # while physical phase and all stiffness values remain unchanged.
    sp_origin, pr_origin = sp_path.epsilon / 2.0, pr_path.epsilon / 2.0
    phase, sp_pair_1, sp_pair_2, ksp, sp_count = mesh_curve(
        sp_path, sun, planet, face_sp, alpha_sp, alpha_nom, young, poisson,
        x_each, x_each, sp_origin)
    _, pr_pair_1, pr_pair_2, kpr, pr_count = mesh_curve(
        pr_path, planet, ring, face_pr, alpha_pr, alpha_nom, young, poisson,
        phase_origin=pr_origin, foundation_second=ring_foundation_compliance)
    values = np.r_[ksp, kpr]
    if not np.all(np.isfinite(values)) or not np.all(values > 0):
        raise RuntimeError("NaN, Inf, or non-positive stiffness encountered.")
    fcarrier = zs / (zs + zr) * sun_rpm / 60
    fmesh = zr * fcarrier
    write_csv("tvms_sp.csv", phase, sp_pair_1, sp_pair_2, ksp, sp_count, "sp", fmesh)
    write_csv("tvms_pr.csv", phase, pr_pair_1, pr_pair_2, kpr, pr_count, "pr", fmesh)
    angles = np.arange(n_planet) * 2 * math.pi / n_planet
    sp_shifts = (zs * angles / (2 * math.pi)) % 1
    pr_shifts = (zr * angles / (2 * math.pi)) % 1
    # The current 21/31/84 three-planet arrangement produces exact zero
    # tooth-index shifts.  No interpolator is used to generate planet data.
    shifted_sp = [ksp.copy() for _ in sp_shifts]
    shifted_pr = [kpr.copy() for _ in pr_shifts]
    with (OUT / "tvms_planets.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["phase_norm", *[f"k_sp_p{i + 1}_N_per_m" for i in range(n_planet)],
                         *[f"k_pr_p{i + 1}_N_per_m" for i in range(n_planet)]])
        for row in zip(phase, *shifted_sp, *shifted_pr):
            writer.writerow([f"{item:.12g}" for item in row])
    configure_matplotlib()
    save_extended_total_figure("Fig2_8a_SP_TVMS.png", phase, ksp, "SP", "#1f4e79")
    save_extended_total_figure("Fig2_8b_PR_TVMS.png", phase, kpr, "PR", "#a64b00")
    save_two_cycle_figure("Fig2_8_TVMS_two_cycles.png", phase, ksp, kpr)
    # Comparison and planet-phase figures retain total stiffness only.
    save_total_figure("Fig2_8_TVMS_compare.png", phase, [(ksp, "Sun--planet", "#1f4e79"), (kpr, "Planet--ring", "#a64b00")])
    colors = ["#1f4e79", "#4f81bd", "#9dc3e6"]
    save_total_figure("Fig2_8_planet_phase.png", phase, [(sig, f"Planet {i + 1} (SP)", colors[i]) for i, sig in enumerate(shifted_sp)])
    endpoint_sp, endpoint_pr = abs(ksp[0] - ksp[-1]) / ksp.mean(), abs(kpr[0] - kpr[-1]) / kpr.mean()
    source_sp, source_pr = source_value("tvms.sunPlanet.contactRatio"), source_value("tvms.ringPlanet.contactRatio")
    phase_ok = np.allclose(sp_shifts, sp_shifts[0], atol=1e-12) and np.allclose(pr_shifts, pr_shifts[0], atol=1e-12)
    sp_double, pr_double = sp_count == 2, pr_count == 2
    sp_sequence, pr_sequence = count_sequence(sp_count), count_sequence(pr_count)
    checks = [("geometry consistency check", abs((zr-zp)*module/2-centre) < 1e-9 and abs((zs+zp)*module/2-centre) > 1e-6),
              ("contact ratio > 1", sp_path.epsilon > 1 and pr_path.epsilon > 1),
              ("stiffness > 0", bool(np.all(values > 0))), ("no NaN/Inf", bool(np.all(np.isfinite(values))),),
              ("double-pair stiffness > typical single-pair stiffness", bool(ksp[sp_double].mean() > ksp[~sp_double].mean() and kpr[pr_double].mean() > kpr[~pr_double].mean())),
              ("complete single -> double -> single sequence", sp_sequence == [1, 2, 1] and pr_sequence == [1, 2, 1]),
              ("cycle start/end mismatch < 1%", endpoint_sp < .01 and endpoint_pr < .01),
              ("phase-shift consistency among planets", phase_ok)]
    summary = ["Healthy potential-energy TVMS", "", "Parameter source: src/pg_parameters.m (all units SI).",
               f"Zs={zs}; Zp={zp}; Zr={zr}; Np={n_planet}; module={module:.9g} m; alpha_nominal={math.degrees(alpha_nom):.6f} deg.",
               f"a_SP=a_PR={centre:.9g} m; alpha_SP={math.degrees(alpha_sp):.6f} deg; alpha_PR={math.degrees(alpha_pr):.6f} deg.",
               f"B_SP={face_sp:.9g} m; B_PR={face_pr:.9g} m; E={young:.9g} Pa; nu={poisson:.9g}; rho={density:.9g} kg/m^3.",
               f"sun RPM={sun_rpm:.9g}; derived mesh frequency={fmesh:.9g} Hz.", "", "Geometry compatibility:",
               f"zero-centre special reference Zs+2Zp={zs+2*zp}; actual Zr={zr}.",
               f"internal nominal centre={(zr-zp)*module/2:.9g} m; external zero-modification centre={(zs+zp)*module/2:.9g} m.",
               f"equivalent SP working-centre tooth-thickness correction: total={x_total:.9g}, symmetric allocation={x_each:.9g} per external gear.",
               "This is a user-authorised modelling assumption, not a measured manufactured profile shift.",
               "SP foundation term: root-section energy continuation (sun/planet fillet/rim remains unavailable).",
               f"PR ring foundation: perforated annular-rim plane-stress FEM; C_f,ring={ring_foundation_compliance:.9e} m/N; k_f,ring={1/ring_foundation_compliance:.9e} N/m.",
               f"Ring FEM: outer radius={rim_args['outer_radius']:.9g} m, 8 through holes of diameter={rim_args['hole_diameter']:.9g} m at radius={rim_args['hole_radius']:.9g} m; outer boundary fixed.",
               f"Ring FEM mesh convergence: coarse C={rim_coarse.compliance_mean:.9e}, fine C={rim_fine.compliance_mean:.9e} m/N, relative change={rim_convergence:.3%}.", "",
               f"display phase origin: SP={sp_origin:.9f}; PR={pr_origin:.9f} (cyclic shift only; no changed stiffness values).",
               f"SP: active-pair sequence={sp_sequence}; displayed double-pair interval={double_interval_text(phase, sp_count)}; remaining intervals are single-pair.",
               f"PR: active-pair sequence={pr_sequence}; displayed double-pair interval={double_interval_text(phase, pr_count)}; remaining intervals are single-pair.",
               f"source-model contact ratios for comparison: SP={source_sp:.9f}, PR={source_pr:.9f}.",
               f"SP ratio difference (computed-source)={sp_path.epsilon-source_sp:+.9f}; PR ratio difference (computed-source)={pr_path.epsilon-source_pr:+.9f}.", "",
               f"k_sp: {stats(ksp)}", f"k_pr: {stats(kpr)}", f"cycle endpoint mismatch: SP={endpoint_sp:.6%}; PR={endpoint_pr:.6%}.",
               "planet mesh-phase shifts: " + ", ".join(f"P{i+1} SP={sp_shifts[i]:.6f}, PR={pr_shifts[i]:.6f}" for i in range(n_planet)), "", "Validation:",
               *[f"{name}: {'PASS' if passed else 'FAIL'}" for name, passed in checks]]
    (OUT / "tvms_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    fem_summary = ["Perforated annular-rim foundation FEM", "", "Plane-stress annulus from ring tooth-root radius to housing-integral outer rim.", "Eight through holes are removed; outer circular rim is fixed.", f"root tooth thickness load patch={ring_root_thickness:.9e} m", f"coarse: nr={rim_coarse.nr}, ntheta={rim_coarse.ntheta}, nodes={rim_coarse.n_nodes}, triangles={rim_coarse.n_triangles}, Cmean={rim_coarse.compliance_mean:.9e} m/N", f"fine: nr={rim_fine.nr}, ntheta={rim_fine.ntheta}, nodes={rim_fine.n_nodes}, triangles={rim_fine.n_triangles}, Cmean={rim_fine.compliance_mean:.9e} m/N", f"fine tooth-location range: Cmin={rim_fine.compliance_min:.9e}, Cmax={rim_fine.compliance_max:.9e} m/N", f"mesh relative change={rim_convergence:.3%}", f"adopted C_f,ring={ring_foundation_compliance:.9e} m/N; adopted k_f,ring={1/ring_foundation_compliance:.9e} N/m", "Limit: exact tooth-root fillet and full 3-D housing geometry are not represented."]
    (OUT / "ring_rim_fem_summary.txt").write_text("\n".join(fem_summary) + "\n", encoding="utf-8")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
