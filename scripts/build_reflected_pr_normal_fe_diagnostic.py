#!/usr/bin/env python3
"""Build a separate FE-asset sensitivity variant for reflected PR load normal.

The current angle registration maps source-local clockwise angle to the
counter-clockwise ANSYS azimuth. Reflection may also reverse the tangential
part of the contact normal. This script tests that alternative without
overwriting or recalibrating the official ANSYS asset.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat


PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "src_py"))
from fe_ring_asset_io import modal_blocks  # noqa: E402

SOURCE = None
ORIGINAL = PROJECT / "results" / "fe_ring_asset" / "ansys_ring_modal_transfer_asset_exact_four_sensor.mat"
OUTPUT = PROJECT / "results" / "fe_planet_q090_23mode_exploratory_20260929" / "fe_asset_reflected_pr_normal_sensitivity.mat"


def run(source: Path, original: Path, output: Path) -> None:
    asset = loadmat(original)
    freq = np.asarray(asset["frequency_Hz"]).reshape(-1)
    ring_xyz, ring_u = modal_blocks(source, 21165, freq.size)
    radius = np.hypot(ring_xyz[:, 0], ring_xyz[:, 1])
    mid_z = 0.5 * (ring_xyz[:, 2].min() + ring_xyz[:, 2].max())
    mid = np.isclose(ring_xyz[:, 2], mid_z, rtol=0, atol=1e-8)
    r_contact = float(asset["contact_radius_m"].item())
    select = mid & np.isclose(radius, r_contact, rtol=0, atol=5e-9)
    xyz = ring_xyz[select]
    u = ring_u[:, select, :]
    theta = np.mod(np.arctan2(xyz[:, 1], xyz[:, 0]), 2.0 * np.pi)
    order = np.argsort(theta)
    theta, xyz, u = theta[order], xyz[order], u[:, order, :]
    if not np.allclose(theta, asset["contact_angle_rad"].reshape(-1), rtol=0, atol=1e-12):
        raise ValueError("Raw ANSYS contact grid does not match the saved official asset.")
    er = np.column_stack((np.cos(theta), np.sin(theta)))
    et = np.column_stack((-np.sin(theta), np.cos(theta)))
    alpha = math.radians(20.0)
    original_normal = -math.sin(alpha) * er + math.cos(alpha) * et
    original_shape = np.einsum("mni,ni->nm", u[:, :, :2], original_normal)
    reference_shape = np.asarray(asset["input_shape_normal"])
    mismatch = np.linalg.norm(original_shape - reference_shape) / np.linalg.norm(reference_shape)
    if mismatch > 1e-10:
        raise ValueError(f"Cannot reproduce official FE shape projection: {mismatch:g}")
    reflected_normal = -math.sin(alpha) * er - math.cos(alpha) * et
    reflected_shape = np.einsum("mni,ni->nm", u[:, :, :2], reflected_normal)
    modified = {key: value for key, value in asset.items() if not key.startswith("__")}
    modified["input_shape_normal"] = reflected_shape
    modified["normal_hypothesis"] = "azimuth_reflection_also_reverses_tangential_PR_normal; diagnostic only"
    output.parent.mkdir(parents=True, exist_ok=True)
    savemat(output, modified, do_compression=True)
    print(f"official projection relative error: {mismatch:.3e}")
    print(f"diagnostic asset: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        required=True,
        help="ANSYS-exported mode_shapes_ring_teeth.csv; raw FE assets stay outside Git.",
    )
    parser.add_argument("--original", type=Path, default=ORIGINAL)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    run(args.source, args.original, args.output)
