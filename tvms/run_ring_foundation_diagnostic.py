"""Read-only sensitivity diagnostic for the current ring foundation term."""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from run_tvms import Gear, OUT, section_width, source_value


def load_rows(name: str) -> list[dict[str, float]]:
    with (OUT / name).open(newline="", encoding="utf-8") as handle:
        raw = list(csv.DictReader(handle))
    rows = []
    for row in raw:
        rows.append({"phase": float(row["phase_norm"]), "count": int(row["active_pair_count"]),
                     "pair1": float(row["k_pair_1_N_per_m"]),
                     "pair2": None if not row["k_pair_2_N_per_m"] else float(row["k_pair_2_N_per_m"]),
                     "total": float(row["k_total_N_per_m"])})
    return rows


def adjusted_pair(k_current: float, c_foundation: float, scale: float | None) -> float:
    """Replace only the existing ring foundation compliance in one pair."""
    c_new = 1.0/k_current - c_foundation if scale is None else 1.0/k_current - c_foundation + c_foundation/scale
    if c_new <= 0:
        raise RuntimeError("Non-positive sensitivity compliance")
    return 1.0/c_new


def main() -> None:
    zs, zp, zr = (int(source_value(x)) for x in ("gear.zs", "gear.zp", "gear.zr"))
    module, alpha, centre = (source_value(x) for x in ("gear.module", "gear.alpha_nominal", "gear.carrierRadius"))
    face, young, poisson = (source_value(x) for x in ("gear.faceWidthRing", "material.gear.E", "material.gear.nu"))
    ring = Gear("ring", zr, module*zr/2, module*zr*math.cos(alpha)/2, source_value("gear.addendumRadiusRing"), source_value("gear.rootRadiusRing"), True)
    planet = Gear("planet", zp, module*zp/2, module*zp*math.cos(alpha)/2, source_value("gear.addendumRadiusPlanet"), source_value("gear.rootRadiusPlanet"))
    alpha_pr = math.acos((ring.pitch-planet.pitch)*math.cos(alpha)/centre)
    shear_modulus = young/(2*(1+poisson))
    thickness = section_width(ring, ring.root, alpha)
    area = face*thickness
    inertia = face*thickness**3/12
    lf_raw = abs(ring.base-ring.root)
    lf = max(lf_raw, 0.05e-3)
    fa, fb = math.sin(alpha_pr), math.cos(alpha_pr)
    cf_bending = fb**2*lf**3/(3*young*inertia)
    cf_axial = fa**2*lf/(young*area)
    cf_shear = 1.2*fb**2*lf/(shear_modulus*area)
    cf = cf_bending+cf_axial+cf_shear
    kf = 1/cf
    pr_rows, sp_rows = load_rows("tvms_pr.csv"), load_rows("tvms_sp.csv")
    sp_mean = float(np.mean([row["total"] for row in sp_rows]))
    single_rows = [row for row in pr_rows if row["count"] == 1]
    typical = min(single_rows, key=lambda row: abs(row["total"]-np.mean([item["total"] for item in single_rows])))
    scales: list[tuple[str, float | None]] = [("0.5x", .5), ("1x", 1.0), ("2x", 2.0), ("5x", 5.0), ("10x", 10.0), ("infinity", None)]
    output = []
    for label, scale in scales:
        totals = []
        for row in pr_rows:
            total = adjusted_pair(row["pair1"], cf, scale)
            if row["pair2"] is not None:
                total += adjusted_pair(row["pair2"], cf, scale)
            totals.append(total)
        pair_typical = adjusted_pair(typical["pair1"], cf, scale)
        output.append({"label": label, "scale": "inf" if scale is None else f"{scale:g}", "kf": "inf" if scale is None else kf*scale,
                       "cf": 0.0 if scale is None else cf/scale, "pair": pair_typical,
                       "mean": float(np.mean(totals)), "ratio": sp_mean/float(np.mean(totals))})
    csv_path = OUT / "ring_foundation_sensitivity.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["case", "foundation_stiffness_scale", "k_f_ring_N_per_m", "C_f_ring_m_per_N", "typical_single_phase", "k_PR_pair_N_per_m", "k_PR_mean_N_per_m", "k_SP_mean_over_k_PR_mean"])
        for row in output:
            writer.writerow([row["label"], row["scale"], row["kf"], f"{row['cf']:.12e}", f"{typical['phase']:.12f}", f"{row['pair']:.12e}", f"{row['mean']:.12e}", f"{row['ratio']:.12f}"])
    project_search = ("Project search found no internal-gear-specific foundation model, rim-thickness model, rim/body flexibility model, or fillet geometry model. "
                      "The only foundation implementation is configs/tvms_config.m: foundationModel='root_section_energy'; the root-crack helper also states that a fillet curve is unavailable.")
    lines = ["RING FOUNDATION STIFFNESS DIAGNOSTIC", "", "This is an offline sensitivity analysis only. The formal TVMS calculation, CSV data and figures were not modified.", "",
             "[Current ring foundation]", "The current implementation uses the same root-section energy continuation form as the external gear, with inward internal-gear coordinate r(x)=r_root-x:",
             "C_f,ring = Fb^2*l_f^3/(3*E*I_root) + Fa^2*l_f/(E*A_root) + 1.2*Fb^2*l_f/(G*A_root)",
             "k_f,ring = 1/C_f,ring; G=E/[2*(1+nu)]; Fa=sin(alpha_w,PR); Fb=cos(alpha_w,PR).", "",
             f"ring root radius = {ring.root:.12g} m", f"ring base radius = {ring.base:.12g} m", f"ring addendum radius (internal tooth tip) = {ring.tip:.12g} m", f"tooth-root thickness = {thickness:.12e} m", f"A_root = B*t_root = {area:.12e} m^2", f"I_root = B*t_root^3/12 = {inertia:.12e} m^4", f"l_f raw = |r_base-r_root| = {lf_raw:.12e} m", f"l_f = max(l_f raw, 0.05 mm) = {lf:.12e} m", f"face width B = {face:.12e} m", f"E = {young:.12e} Pa", f"nu = {poisson:.12f}", f"G = {shear_modulus:.12e} Pa", f"alpha_w,PR = {math.degrees(alpha_pr):.12f} deg", "",
             f"C_f,bending = {cf_bending:.12e} m/N", f"C_f,axial = {cf_axial:.12e} m/N", f"C_f,shear = {cf_shear:.12e} m/N", f"C_f,ring = {cf:.12e} m/N", f"k_f,ring = {kf:.12e} N/m", "",
             "[l_f check]", f"The 0.05 mm floor is NOT active: raw l_f / 0.05 mm = {lf_raw/(.05e-3):.3f}.", "The current l_f is 5.679 mm. Its cubed bending contribution dominates C_f,ring; therefore the large compliance comes from applying the root-section continuation over the full base-to-root radial gap, not from the numerical floor.", "",
             "[Sensitivity: only ring foundation stiffness changes]", f"Typical PR single-pair display phase = {typical['phase']:.12f}", "case | k_PR_pair (N/m) | k_PR_mean (N/m) | k_SP_mean/k_PR_mean"]
    for row in output:
        lines.append(f"{row['label']} | {row['pair']:.12e} | {row['mean']:.12e} | {row['ratio']:.12f}")
    lines += ["", "[Search for a fuller ring-foundation model]", project_search, "", "[Assessment]", "Assessment = likely overestimated compliance.", "Reason: the current term is an external-gear root cantilever/root-section continuation applied inward to the ring tooth. It includes no measured rim thickness, ring-body flexibility, gear-body constraint, tooth-fillet geometry, or internal-gear-specific foundation coefficient. The sensitivity table shows that this single unvalidated term controls the PR mean stiffness level.", "The implementation is algebraically consistent with its stated approximation, but the approximation is not physically well supported for this ring until the rim/fillet/hub geometry or a validated internal-gear foundation model is supplied."]
    text = "\n".join(lines) + "\n"
    (OUT / "ring_foundation_diagnostic.txt").write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
