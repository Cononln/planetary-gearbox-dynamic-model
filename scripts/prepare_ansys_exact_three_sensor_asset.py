"""Build the moving-force/three-radial-sensor modal asset.

The three sensor patches are exported by coordinate rotation inside ANSYS.
No sensor gain, angular window, or response scaling is fitted here.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.io import savemat


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src_py"))

from fe_ring_asset_io import (  # noqa: E402
    SENSOR_IDS, frequencies, modal_blocks, sensor_modes, validate_export_qa,
)


def main():
    parser = argparse.ArgumentParser(
        description="Build a three-sensor ANSYS ring modal transfer asset."
    )
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ring-nodes", type=int, default=20520)
    parser.add_argument("--sensor-nodes", type=int, default=9)
    parser.add_argument("--pitch-radius-m", type=float, default=0.063)
    parser.add_argument("--pressure-angle-deg", type=float, default=20.0)
    args = parser.parse_args()
    if args.ring_nodes <= 0 or args.sensor_nodes <= 0:
        parser.error("--ring-nodes and --sensor-nodes must be positive")
    if args.pitch_radius_m <= 0:
        parser.error("--pitch-radius-m must be positive")

    source = args.source.resolve()
    output = args.output.resolve()
    freq = frequencies(source / "modal_summary.csv")
    export_qa = validate_export_qa(
        source,
        expected_mode_count=freq.size,
        expected_sensor_nodes=args.sensor_nodes,
        check_marker_freshness=True,
    )
    ring_xyz, ring_u = modal_blocks(
        source / "mode_shapes_ring_teeth.csv", args.ring_nodes, freq.size
    )
    sensor_xyz, sensor_u = sensor_modes(
        source / "mode_shapes_sensor.csv", freq.size, args.sensor_nodes
    )

    radius = np.hypot(ring_xyz[:, 0], ring_xyz[:, 1])
    mid_z = 0.5 * (ring_xyz[:, 2].min() + ring_xyz[:, 2].max())
    mid = np.isclose(ring_xyz[:, 2], mid_z, rtol=0.0, atol=1e-8)
    available_radii = np.unique(np.round(radius[mid], 9))
    contact_radius = available_radii[
        np.argmin(np.abs(available_radii - args.pitch_radius_m))
    ]
    contact = mid & np.isclose(radius, contact_radius, rtol=0.0, atol=5e-9)
    if np.count_nonzero(contact) < 84:
        raise RuntimeError("The mid-face pitch-circle FE nodes were not found")

    xyz = ring_xyz[contact]
    u = ring_u[:, contact, :]
    theta = np.mod(np.arctan2(xyz[:, 1], xyz[:, 0]), 2 * np.pi)
    order = np.argsort(theta)
    theta, xyz, u = theta[order], xyz[order], u[:, order, :]
    er = np.column_stack((np.cos(theta), np.sin(theta)))
    et = np.column_stack((-np.sin(theta), np.cos(theta)))
    alpha = math.radians(args.pressure_angle_deg)
    normal_rp = -math.sin(alpha) * er + math.cos(alpha) * et
    input_shape = np.einsum("mni,ni->nm", u[:, :, :2], normal_rp)

    sensor_global = np.empty(3)
    sensor_shape = np.empty((3, freq.size))
    sensor_centres = np.empty((3, 3))
    for index, sensor_id in enumerate(SENSOR_IDS):
        centre = sensor_xyz[sensor_id].mean(axis=0)
        angle = math.atan2(centre[1], centre[0])
        radial = np.array([math.cos(angle), math.sin(angle), 0.0])
        sensor_centres[index] = centre
        sensor_global[index] = np.mod(angle, 2 * np.pi)
        sensor_shape[index] = np.einsum(
            "mni,i->m", sensor_u[sensor_id], radial
        ) / args.sensor_nodes

    # Check each named patch, not just the unordered set of angles: otherwise
    # exchanging SENSOR_120 and SENSOR_240 silently swaps the response labels.
    sensor_geometry_local = np.mod(sensor_global[0] - sensor_global, 2 * np.pi)
    expected_local = np.deg2rad([0.0, 120.0, 240.0])
    geometry_error = np.angle(
        np.exp(1j * (sensor_geometry_local - expected_local))
    )
    if np.max(np.abs(np.rad2deg(geometry_error))) > 2.0:
        raise RuntimeError(
            "Sensor labels do not match clockwise 0/120/240-degree patch centres"
        )
    sensor_local = expected_local

    output.mkdir(parents=True, exist_ok=True)
    asset_path = output / "ansys_ring_modal_transfer_asset_exact_three_sensor.mat"
    savemat(
        asset_path,
        {
            "frequency_Hz": freq.reshape(1, -1),
            "contact_angle_rad": theta.reshape(-1, 1),
            "contact_xyz_m": xyz,
            "input_shape_normal": input_shape,
            "sensor_angle_local_deg": np.rad2deg(sensor_local).reshape(1, -1),
            "sensor_angle_global_deg": np.rad2deg(sensor_global).reshape(1, -1),
            "sensor_shape_radial": sensor_shape,
            "sensor_centre_xyz_m": sensor_centres,
            "mass_normalized": np.array([[1]], dtype=np.uint8),
            "mass_normalization_si_verified": np.array([[0]], dtype=np.uint8),
            "sensor_mode_shapes_exact": np.array([[1]], dtype=np.uint8),
            "contact_radius_m": np.array([[contact_radius]]),
            "midface_z_m": np.array([[mid_z]]),
        },
        do_compression=True,
    )
    report = {
        "source_directory": str(source),
        "source_files_modified": False,
        "sensor_generation": "ANSYS coordinate rotation plus nearest full-mesh node",
        "sensor_observation": (
            f"actual radial displacement averaged over each {args.sensor_nodes}-node patch"
        ),
        "fitted_sensor_gain": False,
        "mode_count": int(freq.size),
        "mass_normalized": True,
        "mass_normalization_si_verified": False,
        "apdl_export_qa": export_qa,
        "sensor_ids": list(SENSOR_IDS),
        "sensor_global_deg": np.rad2deg(sensor_global).tolist(),
        "sensor_local_deg": np.rad2deg(sensor_local).tolist(),
        "sensor_centres_m": sensor_centres.tolist(),
    }
    (output / "ansys_ring_modal_transfer_asset_exact_three_sensor.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    print(asset_path)


if __name__ == "__main__":
    main()
