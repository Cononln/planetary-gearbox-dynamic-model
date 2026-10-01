#!/usr/bin/env python3
"""Decompose the saved 0.90 mm sun-crack FE response by mode and PR input.

This is a diagnostic replay of the *existing* 18-DOF PR-force records. It does
not rerun or change the 18-DOF dynamics, FE asset, or saved sensor response.
For each planet and FE mode, the moving PR-force increment drives the same
exact-ZOH modal equation and static initial condition as the corrected path.
All contributions retain their signs and are evaluated at the peak of the
*total* fault-minus-healthy sensor response in the same event window.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.signal import lfilter, lfilter_zi

from remap_saved_fe_path_q090 import (
    DEFAULT_ASSET,
    DEFAULT_OUT as DEFAULT_CORRECTED,
    DEFAULT_SOURCE_DIR,
    FS_EXPECTED,
    ZETA,
    contact_angles,
    load_asset,
    load_source,
)


DEFAULT_OUT = DEFAULT_CORRECTED.parent / "fe_q090_modal_pr_decomposition_20260929"
SENSOR_INDEX = (0, 1)  # SENSOR_0 and SENSOR_90 in the four-sensor asset
SENSOR_LABEL = (0, 90)
STATE_CENTRES_LOCAL_DEG = np.array([270.0, 180.0, 90.0, 0.0])
EVENT_MARGIN_S = 0.004
RELATIVE_L2_LIMIT = 1e-8
MAX_ABSOLUTE_ERROR_LIMIT_MS2 = 1e-10


def read_reference(path: Path, time: np.ndarray) -> np.ndarray:
    """Read only the existing corrected fault-minus-healthy 0/90 columns."""

    expected = [
        "time_s",
        "q0900um_minus_healthy_sensor_000deg_acceleration_ms2",
        "q0900um_minus_healthy_sensor_090deg_acceleration_ms2",
    ]
    with path.open(newline="", encoding="utf-8") as handle:
        header = next(csv.reader(handle))
    columns = [header.index(name) for name in expected]
    data = np.loadtxt(path, delimiter=",", skiprows=1, usecols=columns)
    if data.shape != (time.size, 3):
        raise ValueError(f"Reference CSV shape {data.shape} does not match {time.size} samples.")
    if not np.allclose(data[:, 0], time, rtol=0.0, atol=1e-14):
        raise ValueError("Reference and source MAT time grids differ.")
    return data[:, 1:]


def read_and_check_events(
    path: Path,
    time: np.ndarray,
    fault_loss: np.ndarray,
    theta_fe: np.ndarray,
    sensor0_fe_deg: float,
) -> list[dict[str, object]]:
    """Re-detect source events and verify the saved corrected event table."""

    with path.open(newline="", encoding="utf-8") as handle:
        saved = [
            row for row in csv.DictReader(handle)
            if row["registration"] == "corrected" and float(row["sensor_angle_local_deg"]) == 0.0
        ]
    total_loss = np.sum(np.maximum(fault_loss, 0.0), axis=1)
    threshold = max(float(total_loss.max()) * 1e-10, np.finfo(float).eps * float(total_loss.max()))
    edge = np.diff(np.r_[False, total_loss > threshold, False].astype(np.int8))
    starts = np.flatnonzero(edge == 1)
    ends = np.flatnonzero(edge == -1) - 1
    complete = (starts > 0) & (ends < time.size - 1)
    starts, ends = starts[complete], ends[complete]
    if len(saved) != len(starts) or len(starts) != 48:
        raise ValueError(f"Expected 48 complete matched events, got {len(starts)} detected and {len(saved)} saved.")

    events: list[dict[str, object]] = []
    for ie, (first, last) in enumerate(zip(starts, ends), start=1):
        source_peak = int(first + np.argmax(total_loss[first:last + 1]))
        planet = int(np.argmax(fault_loss[source_peak, :])) + 1
        fe_angle = float(np.rad2deg(theta_fe[source_peak, planet - 1]) % 360.0)
        # FE polar angles are counter-clockwise; user local angles run
        # clockwise from SENSOR_0, whose exported FE azimuth is ~90 deg.
        local_angle = float((sensor0_fe_deg - fe_angle) % 360.0)
        distance = np.abs((local_angle - STATE_CENTRES_LOCAL_DEG + 180.0) % 360.0 - 180.0)
        state = int(np.argmin(distance)) + 1
        if float(np.min(distance)) > 5.0:
            raise ValueError(f"Event {ie} contact {local_angle:.4f} deg has no expected state.")
        saved_row = saved[ie - 1]
        for key, actual in (
            ("event_number", ie),
            ("active_planet", planet),
            ("source_event_start_time_s", float(time[first])),
            ("source_event_peak_time_s", float(time[source_peak])),
            ("source_event_end_time_s", float(time[last])),
            ("contact_angle_global_deg", fe_angle),
        ):
            if abs(float(saved_row[key]) - float(actual)) > (1e-8 if "angle" in key else 1e-12):
                raise ValueError(f"Saved event {ie} disagrees for {key}: {saved_row[key]} versus {actual}.")
        lo = int(np.searchsorted(time, time[first] - EVENT_MARGIN_S, side="left"))
        hi = int(np.searchsorted(time, time[last] + EVENT_MARGIN_S, side="right"))
        events.append({
            "event_number": ie,
            "state": state,
            "active_planet": planet,
            "first": int(first),
            "last": int(last),
            "source_peak": source_peak,
            "window": slice(lo, hi),
            "local_contact_deg": local_angle,
            "fe_contact_deg": fe_angle,
        })
    counts = np.bincount([int(event["state"]) for event in events], minlength=5)[1:]
    if not np.array_equal(counts, np.full(4, 12)):
        raise ValueError(f"Expected 12 events per local contact state; got {counts.tolist()}.")
    return events


def exact_zoh_acceleration_filter(omega: float, zeta: float, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Transfer Q[k] to eta_ddot[k] for the original exact-ZOH update.

    x[k+1] = Ad*x[k] + Bd*Q[k], eta_ddot[k] = Q[k] + C*x[k].
    The filter initial state is set to the steady response of Q[0], which
    is exactly eta[0]=Q[0]/omega**2 and eta_dot[0]=0.
    """

    a_continuous = np.array([[0.0, 1.0], [-omega * omega, -2.0 * zeta * omega]])
    augmented = np.zeros((3, 3), dtype=float)
    augmented[:2, :2] = a_continuous
    augmented[:2, 2] = (0.0, 1.0)
    matrix = expm(augmented * dt)
    ad, bd = matrix[:2, :2], matrix[:2, 2]
    c = np.array([-omega * omega, -2.0 * zeta * omega])
    trace, det = float(np.trace(ad)), float(np.linalg.det(ad))
    a = np.array([1.0, -trace, det])
    b = np.array([
        1.0,
        -trace + float(c @ bd),
        det + float(c @ (ad - trace * np.eye(2)) @ bd),
    ])
    return b, a


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> Path:
    asset = load_asset(args.asset)
    n = int(round(args.duration_s * args.fs))
    healthy = load_source(args.source_dir / "lpm_response_q0000um.mat", n)
    fault = load_source(args.source_dir / "lpm_response_q0900um.mat", n)
    if not np.isclose(float(healthy["depth_mm"]), 0.0, atol=1e-12):
        raise ValueError("Healthy source MAT depth is not zero.")
    if not np.isclose(float(fault["depth_mm"]), 0.9, atol=1e-12):
        raise ValueError("Fault source MAT depth is not 0.90 mm.")
    time = np.asarray(fault["time"])
    carrier = np.asarray(fault["carrier_angle_rad"])
    if not np.array_equal(time, healthy["time"]) or not np.array_equal(carrier, healthy["carrier_angle_rad"]):
        raise ValueError("Healthy and crack source time/carrier grids differ.")
    if not np.isclose(float(fault["fs"]), args.fs, atol=1e-6, rtol=0.0):
        raise ValueError("Source sampling frequency differs from requested value.")
    delta_force = (
        np.asarray(fault["mesh_force_ring_planet"])
        - np.asarray(healthy["mesh_force_ring_planet"])
    )
    reference = read_reference(args.corrected_dir / "corrected_registration_time_signals.csv", time)
    theta_fe = contact_angles(carrier, asset["sensor_angle_global_deg"], "corrected")
    events = read_and_check_events(
        args.corrected_dir / "corrected_event_metrics.csv",
        time,
        np.asarray(fault["fault_loss"]),
        theta_fe,
        float(asset["sensor_angle_global_deg"][0]),
    )

    n_event, n_sensor, n_source = len(events), 2, 3
    frequency = asset["frequency_hz"]
    n_mode = len(frequency)
    peak_index = np.empty((n_event, n_sensor), dtype=int)
    denominator = np.empty((n_event, n_sensor), dtype=float)
    for ie, event in enumerate(events):
        window = event["window"]
        for isensor in range(n_sensor):
            values = reference[window, isensor]
            peak_index[ie, isensor] = window.start + int(np.argmax(np.abs(values)))
            denominator[ie, isensor] = float(np.dot(values, values))
            if denominator[ie, isensor] <= 0.0:
                raise ValueError(f"Event {ie + 1}, sensor {isensor} has zero fault increment energy.")

    peak_contribution = np.empty((n_event, n_sensor, n_source, n_mode), dtype=float)
    projection_share = np.empty_like(peak_contribution)
    reconstruction = np.zeros((n, n_sensor), dtype=float)
    contact_grid = asset["contact_angle_rad"]
    shape_grid = asset["input_shape_normal"]
    xp = np.concatenate(([contact_grid[-1] - 2.0 * np.pi], contact_grid, [contact_grid[0] + 2.0 * np.pi]))
    dt = float(np.median(np.diff(time)))
    filters = [exact_zoh_acceleration_filter(2.0 * np.pi * f, ZETA, dt) for f in frequency]

    for ip in range(n_source):
        angle = theta_fe[:, ip]
        dforce = delta_force[:, ip]
        for im, (b, a) in enumerate(filters):
            fp = np.concatenate((shape_grid[-1:, im], shape_grid[:, im], shape_grid[:1, im]))
            q = dforce * np.interp(angle, xp, fp)
            eta_ddot, _ = lfilter(b, a, q, zi=lfilter_zi(b, a) * q[0])
            for isensor, original_index in enumerate(SENSOR_INDEX):
                contribution = eta_ddot * asset["sensor_shape_radial"][original_index, im]
                reconstruction[:, isensor] += contribution
                peak_contribution[:, isensor, ip, im] = contribution[peak_index[:, isensor]]
                for ie, event in enumerate(events):
                    window = event["window"]
                    projection_share[ie, isensor, ip, im] = (
                        np.dot(contribution[window], reference[window, isensor])
                        / denominator[ie, isensor]
                    )

    error = reconstruction - reference
    relative_l2 = [float(np.linalg.norm(error[:, i]) / np.linalg.norm(reference[:, i])) for i in range(n_sensor)]
    max_abs = [float(np.max(np.abs(error[:, i]))) for i in range(n_sensor)]
    summed_peak = peak_contribution.sum(axis=(2, 3))
    reference_peak = np.column_stack([reference[peak_index[:, i], i] for i in range(n_sensor)])
    peak_error = float(np.max(np.abs(summed_peak - reference_peak)))
    share_error = float(np.max(np.abs(projection_share.sum(axis=(2, 3)) - 1.0)))
    if max(relative_l2) > RELATIVE_L2_LIMIT or max(max_abs) > MAX_ABSOLUTE_ERROR_LIMIT_MS2:
        raise ValueError(f"Modal/source reconstruction failed: relative L2={relative_l2}, max abs={max_abs}.")
    if peak_error > MAX_ABSOLUTE_ERROR_LIMIT_MS2 or share_error > RELATIVE_L2_LIMIT:
        raise ValueError(f"Peak/share additivity failed: peak error={peak_error}, share error={share_error}.")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    detail_rows: list[dict[str, object]] = []
    mode_rows: list[dict[str, object]] = []
    source_rows: list[dict[str, object]] = []
    state_rows: list[dict[str, object]] = []
    mode_peak = peak_contribution.sum(axis=2)
    mode_share = projection_share.sum(axis=2)
    source_peak = peak_contribution.sum(axis=3)
    source_share = projection_share.sum(axis=3)

    def common(ie: int, isensor: int) -> dict[str, object]:
        event = events[ie]
        index = peak_index[ie, isensor]
        return {
            "event_number": int(event["event_number"]),
            "state": int(event["state"]),
            "local_pr_contact_deg": float(event["local_contact_deg"]),
            "active_planet": int(event["active_planet"]),
            "sensor_local_deg": SENSOR_LABEL[isensor],
            "source_event_start_s": float(time[int(event["first"])]),
            "total_delta_peak_time_s": float(time[index]),
            "total_delta_signed_at_peak_ms2": float(reference[index, isensor]),
            "total_delta_abs_peak_ms2": float(abs(reference[index, isensor])),
        }

    for ie in range(n_event):
        for isensor in range(n_sensor):
            base = common(ie, isensor)
            for ip in range(n_source):
                source_rows.append({
                    **base,
                    "pr_input": ip + 1,
                    "signed_at_total_peak_ms2": float(source_peak[ie, isensor, ip]),
                    "signed_fraction_at_total_peak": float(source_peak[ie, isensor, ip] / reference_peak[ie, isensor]),
                    "projection_share_event_window": float(source_share[ie, isensor, ip]),
                })
                for im in range(n_mode):
                    detail_rows.append({
                        **base,
                        "pr_input": ip + 1,
                        "mode": im + 1,
                        "mode_frequency_hz": float(frequency[im]),
                        "signed_at_total_peak_ms2": float(peak_contribution[ie, isensor, ip, im]),
                        "signed_fraction_at_total_peak": float(peak_contribution[ie, isensor, ip, im] / reference_peak[ie, isensor]),
                        "projection_share_event_window": float(projection_share[ie, isensor, ip, im]),
                    })
            for im in range(n_mode):
                mode_rows.append({
                    **base,
                    "mode": im + 1,
                    "mode_frequency_hz": float(frequency[im]),
                    "signed_at_total_peak_ms2": float(mode_peak[ie, isensor, im]),
                    "signed_fraction_at_total_peak": float(mode_peak[ie, isensor, im] / reference_peak[ie, isensor]),
                    "projection_share_event_window": float(mode_share[ie, isensor, im]),
                })

    for state in range(1, 5):
        members = [ie for ie, event in enumerate(events) if event["state"] == state]
        for isensor in range(n_sensor):
            avg_source_share = source_share[members, isensor, :].mean(axis=0)
            avg_mode_share = mode_share[members, isensor, :].mean(axis=0)
            active_shares = np.array([
                source_share[ie, isensor, int(events[ie]["active_planet"]) - 1]
                for ie in members
            ])
            active_peak_fractions = np.array([
                source_peak[ie, isensor, int(events[ie]["active_planet"]) - 1]
                / reference_peak[ie, isensor]
                for ie in members
            ])
            dominant_modes = np.argsort(-np.abs(avg_mode_share))[:5]
            state_rows.append({
                "state": state,
                "state_centre_local_deg": float(STATE_CENTRES_LOCAL_DEG[state - 1]),
                "sensor_local_deg": SENSOR_LABEL[isensor],
                "event_count": len(members),
                "mean_local_pr_contact_deg": float(np.mean([events[ie]["local_contact_deg"] for ie in members])),
                "mean_total_delta_abs_peak_ms2": float(np.mean(np.abs(reference_peak[members, isensor]))),
                "mean_pr1_projection_share": float(avg_source_share[0]),
                "mean_pr2_projection_share": float(avg_source_share[1]),
                "mean_pr3_projection_share": float(avg_source_share[2]),
                "mean_active_pr_projection_share": float(active_shares.mean()),
                "min_active_pr_projection_share": float(active_shares.min()),
                "max_active_pr_projection_share": float(active_shares.max()),
                "mean_active_pr_signed_fraction_at_total_peak": float(active_peak_fractions.mean()),
                "sum_abs_mean_mode_projection_shares": float(np.sum(np.abs(avg_mode_share))),
                "top5_modes_by_abs_mean_projection_share": ",".join(str(int(im + 1)) for im in dominant_modes),
                "top5_mode_frequency_hz": ",".join(f"{frequency[im]:.6f}" for im in dominant_modes),
                "top5_mode_mean_projection_shares": ",".join(f"{avg_mode_share[im]:.9g}" for im in dominant_modes),
            })

    write_csv(args.out_dir / "event_mode_pr_detail.csv", list(detail_rows[0]), detail_rows)
    write_csv(args.out_dir / "event_mode_sums.csv", list(mode_rows[0]), mode_rows)
    write_csv(args.out_dir / "event_pr_input_sums.csv", list(source_rows[0]), source_rows)
    write_csv(args.out_dir / "state_summary.csv", list(state_rows[0]), state_rows)
    checks = {
        "method": "Existing 18-DOF PR-force increment -> corrected moving FE input-shape -> 60 exact-ZOH modes -> radial 0/90 sensor shapes",
        "source_healthy_mat": str(args.source_dir / "lpm_response_q0000um.mat"),
        "source_fault_mat": str(args.source_dir / "lpm_response_q0900um.mat"),
        "fe_asset_mat": str(args.asset),
        "reference_csv": str(args.corrected_dir / "corrected_registration_time_signals.csv"),
        "sampling_hz": float(fault["fs"]),
        "duration_s": float(args.duration_s),
        "sample_count": n,
        "mode_count": n_mode,
        "source_pr_inputs": n_source,
        "event_count": n_event,
        "events_per_state": [sum(int(event["state"]) == state for event in events) for state in range(1, 5)],
        "static_modal_initial_state": "eta0=Q0/omega^2; eta_dot0=0 for each PR input/mode",
        "modal_damping_ratio": ZETA,
        "sensor_local_deg": list(SENSOR_LABEL),
        "relative_l2_error_0_90": relative_l2,
        "maximum_absolute_error_ms2_0_90": max_abs,
        "maximum_event_peak_additivity_error_ms2": peak_error,
        "maximum_event_projection_share_sum_error": share_error,
        "limits": {
            "relative_l2": RELATIVE_L2_LIMIT,
            "maximum_absolute_error_ms2": MAX_ABSOLUTE_ERROR_LIMIT_MS2,
        },
        "projection_share_definition": "dot(component,total_delta)/dot(total_delta,total_delta) over each original event window; signed and additive",
        "cancellation_index_definition": "sum of absolute per-mode mean projection shares; values above 1 indicate signed cancellation and are not energy percentages",
        "peak_definition": "maximum absolute fault-minus-healthy sensor acceleration within original source event plus/minus 0.004 s",
        "local_contact_angle_definition": "(FE angle of SENSOR_0 - FE PR contact angle) mod 360 degrees, clockwise from SENSOR_0",
        "units_note": "Absolute FE-path acceleration depends on exported ANSYS modal mass normalization; this remains independently unverified.",
    }
    (args.out_dir / "reconstruction_checks.json").write_text(
        json.dumps(checks, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "output": str(args.out_dir),
        "events": n_event,
        "modes": n_mode,
        "relative_l2_0_90": relative_l2,
        "max_abs_error_ms2_0_90": max_abs,
        "peak_additivity_error_ms2": peak_error,
        "projection_share_error": share_error,
    }, ensure_ascii=False, indent=2))
    return args.out_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--asset", type=Path, default=DEFAULT_ASSET)
    parser.add_argument("--corrected-dir", type=Path, default=DEFAULT_CORRECTED)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--duration-s", type=float, default=2.0)
    parser.add_argument("--fs", type=float, default=FS_EXPECTED)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
