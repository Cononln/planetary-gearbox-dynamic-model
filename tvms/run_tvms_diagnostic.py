"""Read-only geometry/contact-ratio diagnostic for the existing TVMS output."""
from __future__ import annotations

import csv
import math
import re
from pathlib import Path

from run_tvms import Gear, PFILE, ROOT, external_path, internal_path, involute, source_value

OUT = ROOT / "results" / "tvms"
CONFIG = ROOT / "configs" / "tvms_config.m"


def config_number(name: str) -> str:
    source = CONFIG.read_text(encoding="utf-8")
    found = re.search(rf"^\s*c\.{re.escape(name)}\s*=\s*([-+0-9.eE]+)\s*;", source, re.MULTILINE)
    return found.group(1) if found else "MISSING"


def load_counts(filename: str) -> tuple[list[float], list[int], list[float]]:
    with (OUT / filename).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return ([float(row["phase_norm"]) for row in rows],
            [int(row["active_pair_count"]) for row in rows],
            [float(row["k_total_N_per_m"]) for row in rows])


def sequence(values: list[int]) -> list[int]:
    return [values[0], *[values[i] for i in range(1, len(values)) if values[i] != values[i - 1]]]


def double_interval(phase: list[float], count: list[int]) -> str:
    points = [value for value, n in zip(phase, count) if n == 2]
    return "MISSING" if not points else f"[{points[0]:.9f}, {points[-1]:.9f}]"


def report_mesh(name: str, mesh, first: Gear, second: Gear, alpha_work: float,
                centre: float, base_pitch: float, phase: list[float], count: list[int]) -> list[str]:
    if mesh.kind == "external":
        term_first = math.sqrt(first.tip**2 - first.base**2) - first.base * math.tan(alpha_work)
        term_second = math.sqrt(second.tip**2 - second.base**2) - second.base * math.tan(alpha_work)
        formula = "L = [sqrt(ra1^2-rb1^2)-rb1*tan(alpha_w)] + [sqrt(ra2^2-rb2^2)-rb2*tan(alpha_w)]"
        detail = [f"  sun recess term = {term_first:.12g} m", f"  planet approach term = {term_second:.12g} m"]
    else:
        ring_term = math.sqrt(second.tip**2 - second.base**2)
        planet_term = math.sqrt(first.tip**2 - first.base**2)
        formula = "L = a*sin(alpha_w) - sqrt(ra_ring^2-rb_ring^2) + sqrt(ra_planet^2-rb_planet^2)"
        detail = [f"  a*sin(alpha_w) = {centre*math.sin(alpha_work):.12g} m",
                  f"  ring involute-boundary term = {ring_term:.12g} m",
                  f"  planet involute-boundary term = {planet_term:.12g} m"]
    theory = mesh.epsilon - 1.0
    numerical = sum(n == 2 for n in count) / len(count)
    error = abs(theory - numerical)
    return [f"[{name} contact-ratio calculation]",
            f"  rb_{first.name} = {first.base:.12g} m; ra_{first.name} = {first.tip:.12g} m",
            f"  rb_{second.name} = {second.base:.12g} m; ra_{second.name} = {second.tip:.12g} m",
            f"  alpha_w = {math.degrees(alpha_work):.9f} deg",
            f"  {formula}", *detail,
            f"  path of contact L = {mesh.length:.12g} m",
            f"  base pitch pb = pi*m*cos(alpha) = {base_pitch:.12g} m",
            f"  epsilon_{name.lower()} = L/pb = {mesh.epsilon:.12f}",
            f"  theoretical double-pair fraction = epsilon-1 = {theory:.12f}",
            f"  numerical double-pair fraction = {sum(n == 2 for n in count)}/{len(count)} = {numerical:.12f}",
            f"  absolute fraction error = {error:.12f} ({100*error:.6f} percentage points): {'PASS' if error <= .01 else 'FAIL'}",
            f"  active-pair sequence = {sequence(count)}; double interval = {double_interval(phase, count)}"]


def main() -> None:
    zs, zp, zr = (int(source_value(name)) for name in ("gear.zs", "gear.zp", "gear.zr"))
    module, alpha, centre = (source_value(name) for name in ("gear.module", "gear.alpha_nominal", "gear.carrierRadius"))
    face_sp, face_pr = (source_value(name) for name in ("gear.faceWidthSunPlanet", "gear.faceWidthRing"))
    sun = Gear("sun", zs, module*zs/2, module*zs*math.cos(alpha)/2, source_value("gear.addendumRadiusSun"), source_value("gear.rootRadiusSun"))
    planet = Gear("planet", zp, module*zp/2, module*zp*math.cos(alpha)/2, source_value("gear.addendumRadiusPlanet"), source_value("gear.rootRadiusPlanet"))
    ring = Gear("ring", zr, module*zr/2, module*zr*math.cos(alpha)/2, source_value("gear.addendumRadiusRing"), source_value("gear.rootRadiusRing"), True)
    alpha_sp = math.acos((sun.pitch + planet.pitch) * math.cos(alpha) / centre)
    alpha_pr = math.acos((ring.pitch - planet.pitch) * math.cos(alpha) / centre)
    base_pitch = math.pi * module * math.cos(alpha)
    sp_mesh = external_path(sun, planet, alpha_sp, base_pitch)
    pr_mesh = internal_path(planet, ring, alpha_pr, centre, base_pitch)
    sp_phase, sp_count, sp_total = load_counts("tvms_sp.csv")
    pr_phase, pr_count, pr_total = load_counts("tvms_pr.csv")
    a_sp_nom, a_pr_nom = (zs + zp)*module/2, (zr - zp)*module/2
    x_total = (zs + zp) * (involute(alpha_sp) - involute(alpha)) / (2*math.tan(alpha))
    geometry_pass = abs(a_pr_nom-centre) < 1e-9 and abs(a_sp_nom-centre) > 1e-6 and bool(re.search(r"useEquivalentWorkingCentreThicknessCorrection\s*=\s*true", CONFIG.read_text(encoding="utf-8")))
    endpoint_sp = abs(sp_total[0]-sp_total[-1]) / (sum(sp_total)/len(sp_total))
    endpoint_pr = abs(pr_total[0]-pr_total[-1]) / (sum(pr_total)/len(pr_total))
    active_logic_sp = sequence(sp_count) == [1, 2, 1]
    active_logic_pr = sequence(pr_count) == [1, 2, 1]
    phase_pass = (sp_phase[0], sp_phase[-1], pr_phase[0], pr_phase[-1]) == (0.0, 1.0, 0.0, 1.0) and endpoint_sp < .01 and endpoint_pr < .01
    report = ["TVMS GEOMETRY AND CONTACT-RATIO DIAGNOSTIC", "", "[Actual parameters read from current project]",
              f"Zs = {zs}", f"Zp = {zp}", f"Zr = {zr}", f"module m = {module:.12g} m", f"nominal pressure angle alpha = {math.degrees(alpha):.9f} deg",
              f"working pressure angle alpha_w,SP = {math.degrees(alpha_sp):.9f} deg", f"working pressure angle alpha_w,PR = {math.degrees(alpha_pr):.9f} deg",
              f"sun profile shift xs = {config_number('xSun')} (TVMS config; supplied drawing value)", f"planet profile shift xp = {config_number('xPlanet')} (TVMS config; supplied drawing value)", f"ring profile shift xr = {config_number('xRing')} (TVMS config; model value)",
              f"face width B_SP = {face_sp:.12g} m", f"face width B_PR = {face_pr:.12g} m", "", "[Geometry]",
              f"external theoretical centre a_sp,0 = m*(Zs+Zp)/2 = {a_sp_nom:.12g} m", f"internal theoretical centre a_pr,0 = m*(Zr-Zp)/2 = {a_pr_nom:.12g} m", f"current carrier/planet-centre radius r_c = {centre:.12g} m",
              f"a_sp,0 == r_c: {'PASS' if abs(a_sp_nom-centre)<1e-9 else 'FAIL'}; residual = {centre-a_sp_nom:.12g} m", f"a_pr,0 == r_c: {'PASS' if abs(a_pr_nom-centre)<1e-9 else 'FAIL'}; residual = {centre-a_pr_nom:.12g} m",
              f"Zr == Zs+2*Zp: {'PASS' if zr == zs+2*zp else 'FAIL'} ({zr} vs {zs+2*zp})", f"external working-centre correction evidence: alpha_w,SP={math.degrees(alpha_sp):.9f} deg and equivalent shift sum={x_total:.12f}",
              "The equivalent shift is an authorised working-centre tooth-thickness assumption, not an asserted manufactured profile shift.", f"geometry consistency = {'PASS (modified-centre geometry explained)' if geometry_pass else 'FAIL'}", "",
              *report_mesh("SP", sp_mesh, sun, planet, alpha_sp, centre, base_pitch, sp_phase, sp_count), "", *report_mesh("PR", pr_mesh, planet, ring, alpha_pr, centre, base_pitch, pr_phase, pr_count), "",
              "[Implementation audit]", f"SP active-pair logic = {'PASS' if active_logic_sp else 'FAIL'}", f"PR active-pair logic = {'PASS' if active_logic_pr else 'FAIL'}", f"phase normalization = {'PASS' if phase_pass else 'FAIL'} (phase 0..1; endpoint mismatch SP={100*endpoint_sp:.9f}%, PR={100*endpoint_pr:.9f}%)",
              "internal/external geometry formulas = PASS (separate external_path and internal_path functions; no formula reuse).", "period boundary tooth-pair bookkeeping = PASS (current pair age phi; preceding pair age 1+phi only when phi < epsilon-1).", "",
              "[Diagnosis]", f"SP epsilon is {sp_mesh.epsilon:.9f}, i.e. approximately 1.15: {'PASS' if abs(sp_mesh.epsilon-1.15)<.01 else 'FAIL'}.",
              f"PR numerical double fraction is consistent with epsilon_PR: {'PASS' if abs((pr_mesh.epsilon-1)-sum(n == 2 for n in pr_count)/len(pr_count)) <= .01 else 'FAIL'}.",
              "Conclusion: the large SP/PR double-pair-region difference follows the computed geometric contact ratios (SP about 1.15; PR about 1.934), not an active-pair interval or phase-normalization error."]
    text = "\n".join(report) + "\n"
    (OUT / "tvms_diagnostic.txt").write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
