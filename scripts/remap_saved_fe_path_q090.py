#!/usr/bin/env python3
"""Re-register a saved 18-DOF PR force record on the FE ring asset.

This is a post-processing diagnostic.  It never calls the 18-DOF solver and
never changes the source force record or the ANSYS asset.  The source record
is a MATLAB v7.3 file, while the small modal asset is a classic MATLAB MAT
file.  The retained source record starts after the 0.2 s discard interval, so
this runner initializes the modal states from the first retained moving-load
sample.  The historical MATLAB run initialized at physical time zero; its
comparison uses a settled interval.

Two registrations are written side by side:

``corrected``
    The agreed convention, with source local zero registered to the exported
    SENSOR_0 global direction and the FE angle increasing counter-clockwise:
    ``theta_fe = sensor0_global - (carrier + planet_phase)``.

``old``
    The historical direct ``theta_fe = carrier + planet_phase`` projection.
    It is retained only so that the saved legacy causal CSV can be checked
    against the Python implementation.

The modal state update is the exact zero-order-hold transition used by
``src/pg_apply_fe_ring_transfer_causal.m``.  Outputs are diagnostic because
the available FE asset's absolute modal normalization remains conditional.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.io import loadmat
from scipy.linalg import expm


DEFAULT_PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_ASSET = (
    DEFAULT_PROJECT
    / "results"
    / "fe_ring_asset"
    / "ansys_ring_modal_transfer_asset_exact_four_sensor.mat"
)
DEFAULT_SOURCE_DIR = (
    DEFAULT_PROJECT / "results" / "four_sensor_lpm_path_comparison_10s_20260923"
)
DEFAULT_OLD_DIR = DEFAULT_PROJECT / "results" / "fe_causal_path_2s_q090_20260926"
DEFAULT_OUT = DEFAULT_PROJECT / "results" / "fe_causal_q090_corrected_registration_20260928"

FS_EXPECTED = 131072.0
ZETA = 0.02
PLANET_PHASE_DEG = np.array([0.0, 120.0, 240.0])


def _scalar(value: np.ndarray) -> float:
    return float(np.asarray(value).reshape(-1)[0])


def load_asset(path: Path) -> dict[str, np.ndarray]:
    """Load the classic MATLAB FE asset and validate its dimensions."""

    data = loadmat(path, squeeze_me=False)
    required = {
        "frequency_Hz",
        "contact_angle_rad",
        "input_shape_normal",
        "sensor_angle_local_deg",
        "sensor_angle_global_deg",
        "sensor_shape_radial",
        "mass_normalized",
    }
    missing = sorted(required.difference(data))
    if missing:
        raise KeyError(f"FE asset is missing fields: {missing}")
    asset = {
        "frequency_hz": np.asarray(data["frequency_Hz"], dtype=float).reshape(-1),
        "contact_angle_rad": np.asarray(data["contact_angle_rad"], dtype=float).reshape(-1),
        "input_shape_normal": np.asarray(data["input_shape_normal"], dtype=float),
        "sensor_angle_local_deg": np.asarray(data["sensor_angle_local_deg"], dtype=float).reshape(-1),
        "sensor_angle_global_deg": np.asarray(data["sensor_angle_global_deg"], dtype=float).reshape(-1),
        "sensor_shape_radial": np.asarray(data["sensor_shape_radial"], dtype=float),
    }
    if not bool(np.asarray(data["mass_normalized"]).reshape(-1)[0]):
        raise ValueError("Only mass-normalized FE assets are supported.")
    n_contact = asset["contact_angle_rad"].size
    n_mode = asset["frequency_hz"].size
    if asset["input_shape_normal"].shape != (n_contact, n_mode):
        raise ValueError("FE contact grid and input-shape dimensions disagree.")
    if asset["sensor_shape_radial"].shape != (
        asset["sensor_angle_local_deg"].size,
        n_mode,
    ):
        raise ValueError("FE sensor-shape dimensions disagree.")
    order = np.argsort(asset["contact_angle_rad"])
    for key in ("contact_angle_rad", "input_shape_normal"):
        asset[key] = asset[key][order] if key == "contact_angle_rad" else asset[key][order, :]
    return asset


def load_source(path: Path, n_samples: int) -> dict[str, np.ndarray | float]:
    """Read one saved ``-struct`` response from a MATLAB v7.3 file."""

    try:
        import h5py
    except ModuleNotFoundError as exc:  # pragma: no cover - environment-specific
        raise RuntimeError(
            "Reading MATLAB v7.3 source records requires h5py; install "
            "the repository requirements first."
        ) from exc

    with h5py.File(path, "r") as handle:
        required = {
            "time",
            "carrierAngleRad",
            "meshForceRingPlanet",
            "faultStiffnessLossSunPlanet",
            "faultStiffnessLossRingPlanet",
        }
        missing = sorted(required.difference(handle.keys()))
        if missing:
            raise KeyError(f"Source response is missing fields: {missing}")
        available = int(handle["time"].shape[-1])
        if n_samples <= 0 or n_samples > available:
            raise ValueError(f"n_samples={n_samples} is outside 1..{available}.")
        # HDF5 stores MATLAB arrays with dimensions reversed.  Restore the
        # MATLAB N-by-3 convention used by the FE transfer functions.
        time = np.asarray(handle["time"][0, :n_samples], dtype=float)
        carrier = np.asarray(handle["carrierAngleRad"][0, :n_samples], dtype=float)
        force = np.asarray(handle["meshForceRingPlanet"][:, :n_samples], dtype=float).T
        sun_loss = np.asarray(handle["faultStiffnessLossSunPlanet"][:, :n_samples], dtype=float).T
        ring_loss = np.asarray(handle["faultStiffnessLossRingPlanet"][:, :n_samples], dtype=float).T
        depth = _scalar(np.asarray(handle["depth"][()])) if "depth" in handle else np.nan
    if force.shape != (n_samples, 3):
        raise ValueError(f"Expected N-by-3 PR force, got {force.shape}.")
    dt = float(np.median(np.diff(time)))
    if not np.allclose(np.diff(time), dt, rtol=1e-10, atol=1e-14):
        raise ValueError("Saved source time is not uniformly sampled.")
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("Saved source time step is invalid.")
    return {
        "time": time,
        "carrier_angle_rad": carrier,
        "mesh_force_ring_planet": force,
        "fault_loss": np.maximum(sun_loss, 0.0) + np.maximum(ring_loss, 0.0),
        "depth_mm": depth,
        "fs": 1.0 / dt,
    }


def periodic_input_shape(
    theta: np.ndarray, contact_angle: np.ndarray, input_shape: np.ndarray
) -> np.ndarray:
    """Interpolate every FE mode on its periodic contact-angle grid."""

    theta = np.mod(np.asarray(theta, dtype=float), 2.0 * np.pi)
    a = np.asarray(contact_angle, dtype=float).reshape(-1)
    v = np.asarray(input_shape, dtype=float)
    xp = np.concatenate(([a[-1] - 2.0 * np.pi], a, [a[0] + 2.0 * np.pi]))
    vp = np.concatenate((v[-1:, :], v, v[:1, :]), axis=0)
    if theta.ndim != 1:
        raise ValueError("Interpolate one planet at a time to bound memory use.")
    out = np.empty((theta.size, v.shape[1]), dtype=float)
    for im in range(v.shape[1]):
        out[:, im] = np.interp(theta, xp, vp[:, im])
    if not np.all(np.isfinite(out)):
        raise ValueError("FE input-shape interpolation returned non-finite values.")
    return out


def contact_angles(
    carrier_angle_rad: np.ndarray,
    sensor_global_deg: np.ndarray,
    registration: str,
) -> np.ndarray:
    """Return N-by-3 FE contact angles for one registration convention."""

    phases = np.deg2rad(PLANET_PHASE_DEG)[None, :]
    carrier = np.asarray(carrier_angle_rad, dtype=float)[:, None]
    if registration == "corrected":
        # Source local zero is the agreed SENSOR_0 direction.  FE atan2
        # angles increase counter-clockwise, hence the reflected subtraction.
        theta = np.deg2rad(float(sensor_global_deg[0])) - carrier - phases
    elif registration == "old":
        theta = carrier + phases
    else:
        raise ValueError(f"Unknown registration: {registration}")
    return np.mod(theta, 2.0 * np.pi)


def modal_zoh(
    time: np.ndarray,
    mesh_force: np.ndarray,
    carrier_angle: np.ndarray,
    asset: dict[str, np.ndarray],
    registration: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Project moving PR force through FE modes with exact ZOH dynamics."""

    frequency = asset["frequency_hz"]
    wn = 2.0 * np.pi * frequency
    zeta = np.full_like(wn, ZETA)
    dt = float(np.median(np.diff(time)))
    theta = contact_angles(carrier_angle, asset["sensor_angle_global_deg"], registration)
    q_modal = np.zeros((time.size, frequency.size), dtype=float)
    for ip in range(PLANET_PHASE_DEG.size):
        psi = periodic_input_shape(
            theta[:, ip], asset["contact_angle_rad"], asset["input_shape_normal"]
        )
        q_modal += mesh_force[:, ip, None] * psi
    n = time.size
    n_mode = wn.size
    n_sensor = asset["sensor_shape_radial"].shape[0]

    ad11 = np.empty(n_mode)
    ad12 = np.empty(n_mode)
    ad21 = np.empty(n_mode)
    ad22 = np.empty(n_mode)
    bd1 = np.empty(n_mode)
    bd2 = np.empty(n_mode)
    for im, (wi, zi) in enumerate(zip(wn, zeta)):
        a = np.array([[0.0, 1.0], [-wi * wi, -2.0 * zi * wi]])
        augmented = np.zeros((3, 3), dtype=float)
        augmented[:2, :2] = a
        augmented[:2, 2] = (0.0, 1.0)
        aug = expm(augmented * dt)
        ad11[im], ad12[im] = aug[0, :2]
        ad21[im], ad22[im] = aug[1, :2]
        bd1[im], bd2[im] = aug[:2, 2]

    sensor_shape = asset["sensor_shape_radial"]
    eta = q_modal[0, :] / (wn * wn)
    eta_dot = np.zeros(n_mode)
    acceleration = np.empty((n, n_sensor), dtype=float)
    eta_ddot = q_modal[0, :] - 2.0 * zeta * wn * eta_dot - wn * wn * eta
    acceleration[0, :] = eta_ddot @ sensor_shape.T
    for k in range(n - 1):
        eta_new = ad11 * eta + ad12 * eta_dot + bd1 * q_modal[k, :]
        eta_dot_new = ad21 * eta + ad22 * eta_dot + bd2 * q_modal[k, :]
        eta, eta_dot = eta_new, eta_dot_new
        eta_ddot = q_modal[k + 1, :] - 2.0 * zeta * wn * eta_dot - wn * wn * eta
        acceleration[k + 1, :] = eta_ddot @ sensor_shape.T
    return acceleration, theta


def event_rows(
    time: np.ndarray,
    fault_loss: np.ndarray,
    delta: np.ndarray,
    theta: np.ndarray,
    sensor_local_deg: np.ndarray,
    sensor_global_deg: np.ndarray,
    registration: str,
) -> list[dict[str, float | int | str]]:
    """Build the source-event and sensor-distance table used by the MATLAB audit."""

    total_loss = np.sum(np.maximum(fault_loss, 0.0), axis=1)
    threshold = max(float(total_loss.max()) * 1e-10, np.finfo(float).eps * float(total_loss.max()))
    active = total_loss > threshold
    edge = np.diff(np.r_[False, active, False].astype(np.int8))
    starts = np.flatnonzero(edge == 1)
    ends = np.flatnonzero(edge == -1) - 1
    complete = (starts > 0) & (ends < time.size - 1)
    starts, ends = starts[complete], ends[complete]
    rows: list[dict[str, float | int | str]] = []
    half_window = 0.004
    global_rad = np.deg2rad(sensor_global_deg)
    for ie, (start, end) in enumerate(zip(starts, ends), start=1):
        ii_event = np.arange(start, end + 1)
        event_idx = int(ii_event[np.argmax(total_loss[ii_event])])
        planet = int(np.argmax(fault_loss[event_idx, :]))
        contact_deg = float(np.rad2deg(theta[event_idx, planet]) % 360.0)
        t0 = float(time[start])
        tpeak = float(time[event_idx])
        tend = float(time[end])
        window = np.flatnonzero((time >= max(0.0, t0 - half_window)) & (time <= min(time[-1], tend + half_window)))
        tt = time[window]
        for sensor_idx, sensor_angle in enumerate(sensor_global_deg):
            x = delta[window, sensor_idx]
            peak_i = int(np.argmax(np.abs(x)))
            peak_abs = float(np.abs(x[peak_i]))
            cumulative = np.cumsum(x * x)
            if cumulative[-1] > 0.0:
                energy_i = int(np.flatnonzero(cumulative >= 0.05 * cumulative[-1])[0])
                energy_lag = float(tt[energy_i] - t0)
            else:
                energy_lag = float("nan")
            crossing = np.flatnonzero((tt >= t0) & (np.abs(x) >= 0.05 * peak_abs))
            first_lag = float(tt[crossing[0]] - t0) if crossing.size else float("nan")
            angular_distance = abs(float(np.rad2deg(np.arctan2(
                np.sin(theta[event_idx, planet] - global_rad[sensor_idx]),
                np.cos(theta[event_idx, planet] - global_rad[sensor_idx]),
            ))))
            rows.append(
                {
                    "registration": registration,
                    "event_number": ie,
                    "active_planet": planet + 1,
                    "contact_angle_global_deg": contact_deg,
                    "sensor_angle_local_deg": float(sensor_local_deg[sensor_idx]),
                    "sensor_angle_global_deg": float(sensor_angle),
                    "contact_sensor_angular_distance_deg": angular_distance,
                    "source_event_start_time_s": t0,
                    "source_event_peak_time_s": tpeak,
                    "source_event_end_time_s": tend,
                    "first_5pct_peak_crossing_lag_s": first_lag,
                    "five_pct_energy_lag_s": energy_lag,
                    "peak_lag_s": float(tt[peak_i] - t0),
                    "increment_peak_abs_ms2": peak_abs,
                }
            )
    return rows


def source_force_event_rows(
    time: np.ndarray,
    fault_loss: np.ndarray,
    mesh_force_delta: np.ndarray,
) -> list[dict[str, float | int]]:
    """Summarize all three PR-force changes for each detected source event.

    ``active_planet`` identifies the planet whose source crack-loss column is
    largest at the event peak.  It does not imply that the other PR forces are
    unchanged: the coupled 18-DOF solution can redistribute load among all
    three meshes.  Keeping these columns beside the FE event table prevents a
    single-contact distance from being overinterpreted as the complete source.
    """

    total_loss = np.sum(np.maximum(fault_loss, 0.0), axis=1)
    threshold = max(float(total_loss.max()) * 1e-10, np.finfo(float).eps * float(total_loss.max()))
    active = total_loss > threshold
    edge = np.diff(np.r_[False, active, False].astype(np.int8))
    starts = np.flatnonzero(edge == 1)
    ends = np.flatnonzero(edge == -1) - 1
    complete = (starts > 0) & (ends < time.size - 1)
    starts, ends = starts[complete], ends[complete]
    rows: list[dict[str, float | int]] = []
    half_window = 0.004
    for ie, (start, end) in enumerate(zip(starts, ends), start=1):
        ii_event = np.arange(start, end + 1)
        event_idx = int(ii_event[np.argmax(total_loss[ii_event])])
        window = np.flatnonzero(
            (time >= max(0.0, time[start] - half_window))
            & (time <= min(time[-1], time[end] + half_window))
        )
        row: dict[str, float | int] = {
            "event_number": ie,
            "active_planet": int(np.argmax(fault_loss[event_idx, :])) + 1,
            "source_event_start_time_s": float(time[start]),
            "source_event_peak_time_s": float(time[event_idx]),
            "source_event_end_time_s": float(time[end]),
        }
        for ip in range(mesh_force_delta.shape[1]):
            x = mesh_force_delta[window, ip]
            row[f"delta_pr{ip + 1}_rms_N"] = float(np.sqrt(np.mean(x * x)))
            row[f"delta_pr{ip + 1}_peak_abs_N"] = float(np.max(np.abs(x)))
        rows.append(row)
    return rows


def write_csv(path: Path, rows: Iterable[dict[str, object]]) -> None:
    rows = list(rows)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_old_csv(path: Path, n_samples: int) -> dict[str, np.ndarray] | None:
    """Read the historical F-drive causal CSV if it is available."""

    if not path.is_file():
        return None
    names = ["000", "090", "120", "240"]
    keys = ["time_s"]
    for case in ("q0000um", "q0900um"):
        for angle in names:
            key = f"{case}_sensor_{angle}deg_acceleration_ms2"
            keys.append(key)
    out: dict[str, np.ndarray] = {key: np.empty(n_samples, dtype=float) for key in keys}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not set(keys).issubset(set(reader.fieldnames or [])):
            return None
        count = 0
        for count, row in enumerate(reader, start=1):
            if count > n_samples:
                break
            for key in keys:
                out[key][count - 1] = float(row[key])
    if count < n_samples:
        return None
    return out


def run(args: argparse.Namespace) -> Path:
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    asset = load_asset(args.asset)
    n_samples = int(round(args.record_duration * args.fs))
    source_paths = {
        "healthy": args.source_dir / "lpm_response_q0000um.mat",
        "fault": args.source_dir / "lpm_response_q0900um.mat",
    }
    source = {case: load_source(path, n_samples) for case, path in source_paths.items()}
    if not np.isclose(float(source["fault"]["fs"]), args.fs, rtol=0.0, atol=1e-6):
        raise ValueError(f"Source fs={source['fault']['fs']} differs from requested fs={args.fs}.")
    time = np.asarray(source["fault"]["time"], dtype=float)
    carrier = np.asarray(source["fault"]["carrier_angle_rad"], dtype=float)
    mesh_healthy = np.asarray(source["healthy"]["mesh_force_ring_planet"], dtype=float)
    mesh_fault = np.asarray(source["fault"]["mesh_force_ring_planet"], dtype=float)
    if not np.allclose(time, source["healthy"]["time"], rtol=0.0, atol=1e-14):
        raise ValueError("Healthy and cracked source time vectors differ.")
    if not np.allclose(carrier, source["healthy"]["carrier_angle_rad"], rtol=0.0, atol=1e-12):
        raise ValueError("Healthy and cracked carrier-angle vectors differ.")

    result: dict[str, dict[str, np.ndarray]] = {}
    for registration in ("corrected", "old"):
        result[registration] = {}
        for case, force in (("healthy", mesh_healthy), ("fault", mesh_fault)):
            y, theta = modal_zoh(time, force, carrier, asset, registration)
            result[registration][case] = y
            if case == "fault":
                result[registration]["theta"] = theta
        result[registration]["delta"] = result[registration]["fault"] - result[registration]["healthy"]

    sensors_local = asset["sensor_angle_local_deg"]
    # Export only the corrected signed response.  The old-angle projection is
    # retained for the numerical validation table below and is not duplicated
    # as a second full-size time-series CSV.
    signal_path = out_dir / "corrected_registration_time_signals.csv"
    with signal_path.open("w", newline="", encoding="utf-8") as handle:
        names = ["time_s"]
        for case in ("healthy", "fault"):
            names.extend(
                f"{case}_sensor_{int(round(angle)):03d}deg_acceleration_ms2"
                for angle in sensors_local
            )
        names.extend(
            f"q0900um_minus_healthy_sensor_{int(round(angle)):03d}deg_acceleration_ms2"
            for angle in sensors_local
        )
        writer = csv.writer(handle)
        writer.writerow(names)
        for k, t in enumerate(time):
            writer.writerow(
                [float(t)]
                + result["corrected"]["healthy"][k, :].tolist()
                + result["corrected"]["fault"][k, :].tolist()
                + result["corrected"]["delta"][k, :].tolist()
            )

    event_rows_all: list[dict[str, object]] = []
    fault_loss = np.asarray(source["fault"]["fault_loss"], dtype=float)
    source_force_rows = source_force_event_rows(
        time, fault_loss, mesh_fault - mesh_healthy
    )
    write_csv(out_dir / "source_pr_force_event_metrics.csv", source_force_rows)
    for registration in ("corrected", "old"):
        rows = event_rows(
            time,
            fault_loss,
            result[registration]["delta"],
            result[registration]["theta"],
            asset["sensor_angle_local_deg"],
            asset["sensor_angle_global_deg"],
            registration,
        )
        write_csv(out_dir / f"{registration}_event_metrics.csv", rows)
        event_rows_all.extend(rows)
    write_csv(out_dir / "event_metrics_both_registrations.csv", event_rows_all)

    summary_rows: list[dict[str, object]] = []
    for registration in ("corrected", "old"):
        for case in ("healthy", "fault"):
            y = result[registration][case]
            for isensor, angle in enumerate(sensors_local):
                summary_rows.append(
                    {
                        "registration": registration,
                        "case": case,
                        "sensor_local_deg": float(angle),
                        "sensor_global_deg": float(asset["sensor_angle_global_deg"][isensor]),
                        "rms_ms2": float(np.sqrt(np.mean(y[:, isensor] ** 2))),
                        "min_ms2": float(np.min(y[:, isensor])),
                        "max_ms2": float(np.max(y[:, isensor])),
                        "peak_to_peak_ms2": float(np.ptp(y[:, isensor])),
                    }
                )
        delta = result[registration]["delta"]
        for isensor, angle in enumerate(sensors_local):
            summary_rows.append(
                {
                    "registration": registration,
                    "case": "fault_minus_healthy",
                    "sensor_local_deg": float(angle),
                    "sensor_global_deg": float(asset["sensor_angle_global_deg"][isensor]),
                    "rms_ms2": float(np.sqrt(np.mean(delta[:, isensor] ** 2))),
                    "min_ms2": float(np.min(delta[:, isensor])),
                    "max_ms2": float(np.max(delta[:, isensor])),
                    "peak_to_peak_ms2": float(np.ptp(delta[:, isensor])),
                }
            )
    write_csv(out_dir / "response_summary.csv", summary_rows)

    old_csv = load_old_csv(args.old_csv, n_samples)
    validation: list[dict[str, object]] = []
    if old_csv is not None:
        if not np.allclose(time, old_csv["time_s"], rtol=0.0, atol=1e-12):
            raise ValueError("Saved 18-DOF source and legacy causal CSV use different time grids.")
        support_masks = {
            "full_record": np.ones(n_samples, dtype=bool),
            "post_first_event_transient": old_csv["time_s"] >= args.validation_start_s,
        }
        for case in ("healthy", "fault"):
            for isensor, angle in enumerate(sensors_local):
                expected = old_csv[f"q{'0000' if case == 'healthy' else '0900'}um_sensor_{int(round(angle)):03d}deg_acceleration_ms2"]
                predicted = result["old"][case][:, isensor]
                for support, mask in support_masks.items():
                    err = predicted[mask] - expected[mask]
                    denom = max(float(np.linalg.norm(expected[mask])), np.finfo(float).tiny)
                    validation.append(
                        {
                            "support": support,
                            "support_start_s": float(old_csv["time_s"][mask][0]),
                            "case": case,
                            "sensor_local_deg": float(angle),
                            "max_abs_error_ms2": float(np.max(np.abs(err))),
                            "rms_error_ms2": float(np.sqrt(np.mean(err * err))),
                            "relative_l2_error": float(np.linalg.norm(err) / denom),
                            "corrcoef": float(np.corrcoef(predicted[mask], expected[mask])[0, 1]),
                        }
                    )
    write_csv(out_dir / "old_angle_validation_vs_legacy_causal.csv", validation)
    validation_post = [
        row for row in validation if row["support"] == "post_first_event_transient"
    ]
    validation_max_rel_l2 = (
        max(float(row["relative_l2_error"]) for row in validation_post)
        if validation_post
        else float("nan")
    )
    validation_pass = bool(validation_post) and validation_max_rel_l2 <= args.validation_relative_l2_tolerance

    # A compact visual audit is useful for checking the event ordering while
    # retaining signed amplitudes.  Plotting is optional so the numerical
    # result remains usable on a headless installation without Matplotlib.
    try:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 1, figsize=(12, 7), constrained_layout=True)
        for isensor, angle in enumerate(sensors_local):
            axes[0].plot(time, result["corrected"]["delta"][:, isensor] * 1e3, lw=0.6, label=f"SENSOR_{int(round(angle))}")
        axes[0].set_ylabel("fault increment / 10^-3 m s^-2")
        axes[0].set_title("Corrected FE registration: signed fault increment")
        axes[0].grid(True, alpha=0.25)
        axes[0].legend(ncol=4, fontsize=8)
        corrected_events = [r for r in event_rows_all if r["registration"] == "corrected"]
        for sensor_angle in sensors_local:
            rows = [r for r in corrected_events if r["sensor_angle_local_deg"] == float(sensor_angle)]
            axes[1].plot(
                [r["event_number"] for r in rows],
                [r["increment_peak_abs_ms2"] * 1e3 for r in rows],
                ".-",
                ms=3,
                lw=0.7,
                label=f"SENSOR_{int(round(sensor_angle))}",
            )
        axes[1].set_xlabel("source event number")
        axes[1].set_ylabel("event peak / 10^-3 m s^-2")
        axes[1].set_title("Corrected event peak magnitude (no channel scaling)")
        axes[1].grid(True, alpha=0.25)
        axes[1].legend(ncol=4, fontsize=8)
        fig.savefig(out_dir / "corrected_registration_event_diagnostic.png", dpi=180)
        plt.close(fig)
    except Exception as exc:  # pragma: no cover - optional presentation output
        (out_dir / "plot_note.txt").write_text(f"Plot skipped: {exc}\n", encoding="utf-8")

    metadata = {
        "source_healthy": str(source_paths["healthy"]),
        "source_fault": str(source_paths["fault"]),
        "asset": str(args.asset),
        "legacy_causal_csv": str(args.old_csv),
        "n_samples": n_samples,
        "fs_hz": float(source["fault"]["fs"]),
        "time_start_s": float(time[0]),
        "time_end_s": float(time[-1]),
        "carrier_angle_start_rad": float(carrier[0]),
        "carrier_angle_end_rad": float(carrier[-1]),
        "modal_damping_ratio": ZETA,
        "registration_corrected": "theta_fe = sensor0_global - (carrier + planet_phase), modulo 2pi",
        "registration_old": "theta_fe = carrier + planet_phase, modulo 2pi",
        "planet_phase_deg": PLANET_PHASE_DEG.tolist(),
        "initial_modal_state": "static eta[0]=Q[0]/omega^2, eta_dot[0]=0 at first retained sample",
        "source_record_note": "Saved per-depth MAT starts after the original 0.2 s discard; no 18DOF solve was rerun.",
        "legacy_validation_start_s": args.validation_start_s,
        "legacy_validation_relative_l2_tolerance": args.validation_relative_l2_tolerance,
        "legacy_validation_max_relative_l2": validation_max_rel_l2,
        "legacy_validation_pass": validation_pass,
        "legacy_validation_note": "Old-angle projection is compared with the historical MATLAB causal CSV after retained-sample modal settling.",
        "asset_sensor_local_deg": asset["sensor_angle_local_deg"].tolist(),
        "asset_sensor_global_deg": asset["sensor_angle_global_deg"].tolist(),
        "asset_mode_count": int(asset["frequency_hz"].size),
        "asset_frequency_hz_min": float(np.min(asset["frequency_hz"])),
        "asset_frequency_hz_max": float(np.max(asset["frequency_hz"])),
        "healthy_depth_mm": float(source["healthy"]["depth_mm"]),
        "fault_depth_mm": float(source["fault"]["depth_mm"]),
        "source_pr_force_event_metrics": str(out_dir / "source_pr_force_event_metrics.csv"),
        "source_pr_force_note": (
            "All three PR force increments are retained; active_planet identifies "
            "the source crack-loss maximum only."
        ),
    }
    (out_dir / "README.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    (out_dir / "README.txt").write_text(
        "Saved-source FE path re-registration diagnostic.\n"
        "Corrected four-sensor signed responses and per-event metrics are exported.\n"
        "The old-angle projection is used only to validate the saved-source calculation\n"
        "against the historical MATLAB causal response CSV.\n"
        f"Legacy old-angle post-transient validation: {'PASS' if validation_pass else 'FAIL'}; "
        f"max relative L2={validation_max_rel_l2:.9g}; tolerance={args.validation_relative_l2_tolerance:.9g}.\n"
        "The source 18-DOF record and ANSYS asset are read-only inputs. See README.json for provenance.\n",
        encoding="utf-8",
    )
    if not validation_pass:
        raise RuntimeError(
            "Old-angle projection did not match the historical MATLAB causal CSV "
            f"after {args.validation_start_s:g} s (max relative L2 "
            f"{validation_max_rel_l2:.6g}); inspect validation CSV before using corrected results."
        )
    return out_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asset", type=Path, default=DEFAULT_ASSET)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--old-csv", type=Path, default=DEFAULT_OLD_DIR / "causal_fe_path_time_signals.csv")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--fs", type=float, default=FS_EXPECTED)
    parser.add_argument("--record-duration", type=float, default=2.0)
    parser.add_argument(
        "--validation-start-s",
        type=float,
        default=0.02,
        help="Beginning of the old-angle comparison after retained-record modal settling.",
    )
    parser.add_argument(
        "--validation-relative-l2-tolerance",
        type=float,
        default=1e-3,
        help="Maximum allowed old-angle relative L2 error after settling.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    output = run(parse_args())
    print(f"Wrote saved-source FE registration diagnostic to {output}")
