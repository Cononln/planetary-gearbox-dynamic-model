"""Read-only stiffness-component diagnostic for the existing potential-energy TVMS."""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from run_tvms import (Gear, OUT, contact_radii, external_path, internal_path,
                      involute, section_width, source_value)


def terms(gear: Gear, r_contact: float, face_width: float, alpha_work: float,
          alpha_nom: float, young: float, poisson: float, x_equiv: float = 0.0) -> dict[str, float]:
    """Exact compliance terms used by run_tvms.tooth_compliance, unmodified."""
    shear_modulus = young / (2.0 * (1.0 + poisson))
    h = (gear.root-r_contact) if gear.internal else (r_contact-gear.root)
    h = max(h, 0.05e-3)
    fa, fb = math.sin(alpha_work), math.cos(alpha_work)
    def section_at(x: float) -> tuple[float, float]:
        radius = gear.root-x if gear.internal else gear.root+x
        thickness = section_width(gear, radius, alpha_nom, x_equiv)
        return face_width*thickness, face_width*thickness**3/12.0
    cb = quad(lambda x: (fb*(h-x)-fa*h)**2/(young*section_at(x)[1]), 0, h, epsabs=1e-18, epsrel=1e-8, limit=100)[0]
    ca = quad(lambda x: fa**2/(young*section_at(x)[0]), 0, h, epsabs=1e-18, epsrel=1e-8, limit=100)[0]
    cs = quad(lambda x: 1.2*fb**2/(shear_modulus*section_at(x)[0]), 0, h, epsabs=1e-18, epsrel=1e-8, limit=100)[0]
    area0, inertia0 = section_at(0.0)
    lf = max(abs(gear.base-gear.root), 0.05e-3)
    cf = fb**2*lf**3/(3*young*inertia0) + fa**2*lf/(young*area0) + 1.2*fb**2*lf/(shear_modulus*area0)
    return {"bending": cb, "axial": ca, "shear": cs, "foundation": cf}


def hertz_compliance(young: float, poisson: float, face_width: float) -> float:
    equivalent_modulus = 1.0/(2.0*(1.0-poisson**2)/young)
    return 1.0/(math.pi*face_width*equivalent_modulus)


def load_mesh_csv(name: str) -> list[dict[str, float]]:
    with (OUT / name).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return [{"phase": float(row["phase_norm"]), "count": int(row["active_pair_count"]), "total": float(row["k_total_N_per_m"])} for row in rows]


def representative(rows: list[dict[str, float]], count: int) -> dict[str, float]:
    subset = [row for row in rows if row["count"] == count]
    mean = sum(row["total"] for row in subset)/len(subset)
    return min(subset, key=lambda row: abs(row["total"]-mean))


def pair_at(mesh, first: Gear, second: Gear, alpha_work: float, alpha_nom: float,
            young: float, poisson: float, face_width: float, phase: float,
            phase_origin: float, x_first: float = 0.0, x_second: float = 0.0) -> list[dict[str, object]]:
    physical = (phase+phase_origin) % 1.0
    ages = [physical]
    if physical < mesh.epsilon-1.0:
        ages.append(physical+1.0)
    ch = hertz_compliance(young, poisson, face_width)
    result = []
    for age in ages:
        r1, r2 = contact_radii(mesh, first, second, age/mesh.epsilon, alpha_work)
        first_terms, second_terms = terms(first, r1, face_width, alpha_work, alpha_nom, young, poisson, x_first), terms(second, r2, face_width, alpha_work, alpha_nom, young, poisson, x_second)
        c_total = sum(first_terms.values()) + sum(second_terms.values()) + ch
        result.append({"age": age, "r1": r1, "r2": r2, "first": first_terms, "second": second_terms, "hertz": ch, "total_compliance": c_total, "k_pair": 1.0/c_total})
    return result


def stiffness_line(label: str, compliance: float) -> str:
    return f"  {label}: C=1/k={compliance:.12e} m/N; k={1/compliance:.12e} N/m"


def describe_pair(pair: dict[str, object], first: Gear, second: Gear, index: int) -> list[str]:
    first_terms = pair["first"]
    second_terms = pair["second"]
    lines = [f"Pair {index}: contact age={pair['age']:.12f}; r_{first.name}={pair['r1']:.12g} m; r_{second.name}={pair['r2']:.12g} m"]
    for gear, values in ((first, first_terms), (second, second_terms)):
        lines.append(f" {gear.name} tooth:")
        lines += [stiffness_line(f"k_b,{gear.name}", values["bending"]), stiffness_line(f"k_a,{gear.name}", values["axial"]), stiffness_line(f"k_s,{gear.name}", values["shear"]), stiffness_line(f"k_f,{gear.name}", values["foundation"])]
    lines += [stiffness_line("k_h (shared Hertz contact)", pair["hertz"]), f"  total pair compliance = {pair['total_compliance']:.12e} m/N", f"  k_pair,{index} = {pair['k_pair']:.12e} N/m"]
    return lines


def aggregate_categories(pair: dict[str, object]) -> dict[str, float]:
    first, second = pair["first"], pair["second"]
    return {"bending": first["bending"]+second["bending"], "axial": first["axial"]+second["axial"], "shear": first["shear"]+second["shear"], "foundation": first["foundation"]+second["foundation"], "Hertz": pair["hertz"]}


def main() -> None:
    zs, zp, zr = (int(source_value(x)) for x in ("gear.zs", "gear.zp", "gear.zr"))
    module, alpha, centre = (source_value(x) for x in ("gear.module", "gear.alpha_nominal", "gear.carrierRadius"))
    face_sp, face_pr, young, poisson = (source_value(x) for x in ("gear.faceWidthSunPlanet", "gear.faceWidthRing", "material.gear.E", "material.gear.nu"))
    sun = Gear("sun", zs, module*zs/2, module*zs*math.cos(alpha)/2, source_value("gear.addendumRadiusSun"), source_value("gear.rootRadiusSun"))
    planet = Gear("planet", zp, module*zp/2, module*zp*math.cos(alpha)/2, source_value("gear.addendumRadiusPlanet"), source_value("gear.rootRadiusPlanet"))
    ring = Gear("ring", zr, module*zr/2, module*zr*math.cos(alpha)/2, source_value("gear.addendumRadiusRing"), source_value("gear.rootRadiusRing"), True)
    alpha_sp = math.acos((sun.pitch+planet.pitch)*math.cos(alpha)/centre)
    alpha_pr = math.acos((ring.pitch-planet.pitch)*math.cos(alpha)/centre)
    base_pitch = math.pi*module*math.cos(alpha)
    sp_mesh, pr_mesh = external_path(sun, planet, alpha_sp, base_pitch), internal_path(planet, ring, alpha_pr, centre, base_pitch)
    x_total = (zs+zp)*(involute(alpha_sp)-involute(alpha))/(2*math.tan(alpha))
    sp_rows, pr_rows = load_mesh_csv("tvms_sp.csv"), load_mesh_csv("tvms_pr.csv")
    sp_single, sp_double = representative(sp_rows, 1), representative(sp_rows, 2)
    pr_single, pr_double = representative(pr_rows, 1), representative(pr_rows, 2)
    sp_single_pairs = pair_at(sp_mesh, sun, planet, alpha_sp, alpha, young, poisson, face_sp, sp_single["phase"], sp_mesh.epsilon/2, x_total/2, x_total/2)
    sp_double_pairs = pair_at(sp_mesh, sun, planet, alpha_sp, alpha, young, poisson, face_sp, sp_double["phase"], sp_mesh.epsilon/2, x_total/2, x_total/2)
    pr_single_pairs = pair_at(pr_mesh, planet, ring, alpha_pr, alpha, young, poisson, face_pr, pr_single["phase"], pr_mesh.epsilon/2)
    pr_double_pairs = pair_at(pr_mesh, planet, ring, alpha_pr, alpha, young, poisson, face_pr, pr_double["phase"], pr_mesh.epsilon/2)
    sp_categories, pr_categories = aggregate_categories(sp_single_pairs[0]), aggregate_categories(pr_single_pairs[0])
    category_delta = {name: pr_categories[name]-sp_categories[name] for name in sp_categories}
    major = max(category_delta, key=category_delta.get)
    def stats(rows):
        all_k = np.array([row["total"] for row in rows]); single = np.array([row["total"] for row in rows if row["count"] == 1]); double = np.array([row["total"] for row in rows if row["count"] == 2])
        return all_k.min(), all_k.max(), all_k.mean(), single.mean(), double.mean()
    sp_stats, pr_stats = stats(sp_rows), stats(pr_rows)
    lines = ["TVMS STIFFNESS COMPONENT DIAGNOSTIC", "", "No TVMS formula, load-sharing rule, parameter, or curve was changed.", "",
             "[Current foundation-compliance implementation]", "External tooth: r(x)=r_root+x; l_f=max(|r_base-r_root|, 0.05 mm).", "Internal ring tooth: r(x)=r_root-x; l_f=max(|r_base-r_root|, 0.05 mm).", "Both use C_f=Fb^2*l_f^3/(3*E*I_root)+Fa^2*l_f/(E*A_root)+1.2*Fb^2*l_f/(G*A_root).", "Thus the internal ring foundation compliance is included with inward radial integration and its own base-to-root distance.", "It remains a root-section continuation approximation, not a measured ring-rim/hub flexibility model.", "",
             "[SP typical single-pair phase]", f"display phase={sp_single['phase']:.12f}; total TVMS={sp_single['total']:.12e} N/m", *describe_pair(sp_single_pairs[0], sun, planet, 1), "",
             "[SP typical double-pair phase]", f"display phase={sp_double['phase']:.12f}; total TVMS CSV={sp_double['total']:.12e} N/m", *describe_pair(sp_double_pairs[0], sun, planet, 1), *describe_pair(sp_double_pairs[1], sun, planet, 2), f"  k_total = k_pair,1 + k_pair,2 = {sum(pair['k_pair'] for pair in sp_double_pairs):.12e} N/m", "",
             "[PR typical single-pair phase]", f"display phase={pr_single['phase']:.12f}; total TVMS={pr_single['total']:.12e} N/m", *describe_pair(pr_single_pairs[0], planet, ring, 1), "",
             "[PR typical double-pair phase]", f"display phase={pr_double['phase']:.12f}; total TVMS CSV={pr_double['total']:.12e} N/m", *describe_pair(pr_double_pairs[0], planet, ring, 1), *describe_pair(pr_double_pairs[1], planet, ring, 2), f"  k_total = k_pair,1 + k_pair,2 = {sum(pair['k_pair'] for pair in pr_double_pairs):.12e} N/m", "",
             "[TVMS amplitude summary]", f"SP minimum stiffness = {sp_stats[0]:.12e} N/m", f"SP maximum stiffness = {sp_stats[1]:.12e} N/m", f"SP mean stiffness = {sp_stats[2]:.12e} N/m", f"SP single-pair typical stiffness (mean over active_count=1) = {sp_stats[3]:.12e} N/m", f"SP double-pair typical stiffness (mean over active_count=2) = {sp_stats[4]:.12e} N/m", f"PR minimum stiffness = {pr_stats[0]:.12e} N/m", f"PR maximum stiffness = {pr_stats[1]:.12e} N/m", f"PR mean stiffness = {pr_stats[2]:.12e} N/m", f"PR single-pair typical stiffness (mean over active_count=1) = {pr_stats[3]:.12e} N/m", f"PR double-pair typical stiffness (mean over active_count=2) = {pr_stats[4]:.12e} N/m", f"k_SP_mean / k_PR_mean = {sp_stats[2]/pr_stats[2]:.12f}", "",
             "[Typical single-pair compliance composition]", "SP:"]
    for name, value in sp_categories.items():
        lines.append(f"  {name}: {value:.12e} m/N ({100*value/sum(sp_categories.values()):.3f}%)")
    lines.append("PR:")
    for name, value in pr_categories.items():
        lines.append(f"  {name}: {value:.12e} m/N ({100*value/sum(pr_categories.values()):.3f}%; PR-SP delta={category_delta[name]:.12e} m/N)")
    lines += ["", "[Diagnosis]", f"Largest positive PR-minus-SP compliance category at the representative single-pair state: {major}.", "The dominant reason PR is less stiff in the current model is the ring-tooth foundation/root-section compliance; bending flexibility is secondary.", "Hertz compliance is identical for SP and PR because E, nu and face width are identical; it cannot explain their mean-level difference.", "The internal-mesh geometry affects the ring contact/root coordinate and makes its foundation continuation much larger. This is implemented consistently, but its physical magnitude remains sensitive to unavailable fillet/rim/hub geometry."]
    text = "\n".join(lines) + "\n"
    (OUT / "tvms_stiffness_component_diagnostic.txt").write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
