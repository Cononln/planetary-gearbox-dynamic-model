#!/usr/bin/env python3
"""Exploratory four-sensor FE-path audit for saved 18-DOF planet-crack data.

The source pair is the Chapter-2 potential-energy-TVMS healthy/planet-q0.90
run, not the earlier sun-crack LPM run. The 51.2 kHz source can resolve only
the FE modes below its 25.6 kHz Nyquist frequency. No source dynamics, FE
shapes, waveforms, or gains are altered. SP and PR fault events are kept
separate, while all three PR forces drive one common four-sensor path model.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.signal import lfilter, lfilter_zi

from decompose_fe_q090_modal_pr_contributions import exact_zoh_acceleration_filter
from remap_saved_fe_path_q090 import contact_angles, load_asset, modal_zoh


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "results" / "thesis_ch2_fault_dynamics_20260912_094050"
ASSET = PROJECT / "results" / "fe_ring_asset" / "ansys_ring_modal_transfer_asset_exact_four_sensor.mat"
OUTPUT = PROJECT / "results" / "fe_planet_q090_23mode_exploratory_20260929"
ZETA = 0.02
SETTLED_START_S = 0.2
EVENT_TAIL_S = 0.004


def read_state(path: Path) -> dict[str, np.ndarray | float]:
    raw = loadmat(path)
    fields = {
        "time": "time_s",
        "phase": "mesh_phase_cycles",
        # Match pg_apply_fe_ring_transfer_causal, which takes the source
        # solver's signed meshForceRingPlanet without a sign flip.
        "force": "mesh_force_pr_signed_N",
        "sp_loss": "fault_stiffness_loss_sp_N_per_m",
        "pr_loss": "fault_stiffness_loss_pr_N_per_m",
    }
    state: dict[str, np.ndarray | float] = {}
    for key, matlab_key in fields.items():
        if matlab_key not in raw:
            raise KeyError(f"{path} lacks {matlab_key}")
        state[key] = np.asarray(raw[matlab_key], dtype=float).reshape(-1) if key in ("time", "phase") else np.asarray(raw[matlab_key], dtype=float)
    state["fs"] = float(raw["fs_Hz"].item())
    state["dof"] = int(raw["model_DOF"].item())
    if not np.allclose(raw["mesh_force_pr_compression_N"], -raw["mesh_force_pr_signed_N"], rtol=0, atol=1e-10):
        raise ValueError("Compression and signed PR-force conventions disagree.")
    if np.asarray(state["force"]).shape != (np.asarray(state["time"]).size, 3):
        raise ValueError(f"Invalid three-PR-force shape in {path}")
    return state


def event_segments(time: np.ndarray, loss: np.ndarray, kind: str) -> list[dict[str, object]]:
    """Detect complete target-planet tooth encounters after the transient."""
    if np.any(loss[:, 1:] > 1e-6):
        raise ValueError(f"{kind} loss acts outside the selected P1 fault branch.")
    active = loss[:, 0] > max(float(np.max(loss[:, 0])) * 1e-10, 1e-6)
    edges = np.diff(np.r_[False, active, False].astype(np.int8))
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1) - 1
    result: list[dict[str, object]] = []
    for first, last in zip(starts, ends):
        if time[first] < SETTLED_START_S or last >= len(time) - 1:
            continue
        if time[last] + EVENT_TAIL_S > time[-1]:
            continue
        source_peak = first + int(np.argmax(loss[first : last + 1, 0]))
        window_end = int(np.searchsorted(time, time[last] + EVENT_TAIL_S, side="right"))
        result.append({
            "kind": kind,
            "first": int(first),
            "last": int(last),
            "peak": int(source_peak),
            "window": slice(int(first), window_end),
        })
    return result


def save_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty diagnostic table: {path}")
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> dict[str, object]:
    healthy_path = args.source / "case_healthy_state.mat"
    fault_path = args.source / "case_planet_root_crack_qc0p900mm_state.mat"
    healthy, fault = read_state(healthy_path), read_state(fault_path)
    time = np.asarray(fault["time"])
    phase = np.asarray(fault["phase"])
    force_h = np.asarray(healthy["force"])
    force_f = np.asarray(fault["force"])
    sp_loss = np.asarray(fault["sp_loss"])
    pr_loss = np.asarray(fault["pr_loss"])
    fs = float(fault["fs"])
    if healthy["dof"] != 18 or fault["dof"] != 18:
        raise ValueError("Source cases are not both 18 DOF.")
    if not np.array_equal(time, healthy["time"]) or not np.array_equal(phase, healthy["phase"]):
        raise ValueError("Healthy and fault sources do not share the time/phase grid.")
    if not np.isclose(fs, 51200.0, rtol=0, atol=1e-8) or not np.isclose(float(healthy["fs"]), fs):
        raise ValueError("Source sampling rate changed; re-audit mode selection.")
    if not np.allclose(np.diff(time), 1.0 / fs, rtol=0, atol=1e-13):
        raise ValueError("Nonuniform source time grid.")
    if np.max(np.abs(np.asarray(healthy["sp_loss"]))) > 1e-8 or np.max(np.abs(np.asarray(healthy["pr_loss"]))) > 1e-8:
        raise ValueError("Healthy source contains a nonzero fault stiffness loss.")

    asset = load_asset(args.asset)
    all_frequency = np.asarray(asset["frequency_hz"])
    keep = all_frequency < fs / 2.0
    if not np.any(keep):
        raise ValueError("No FE mode is below source Nyquist frequency.")
    frequency = all_frequency[keep]
    shape_in = np.asarray(asset["input_shape_normal"])[:, keep]
    shape_sensor = np.asarray(asset["sensor_shape_radial"])[:, keep]
    angle_grid = np.asarray(asset["contact_angle_rad"])
    if frequency.size != 23 or all_frequency.size != 60:
        raise ValueError("Expected the reviewed 23-of-60 FE mode selection.")
    sensor_local = np.asarray(asset["sensor_angle_local_deg"])
    sensor_global = np.asarray(asset["sensor_angle_global_deg"])
    if not np.allclose(sensor_local, [0, 90, 120, 240]):
        raise ValueError("Unexpected sensor ordering in FE asset.")

    # The saved source uses mesh_counter = Zr*carrier_angle/(2*pi), Zr=84.
    carrier = 2.0 * np.pi * phase / 84.0
    theta_fe = contact_angles(carrier, sensor_global, "corrected")
    xp = np.r_[angle_grid[-1] - 2.0 * np.pi, angle_grid, angle_grid[0] + 2.0 * np.pi]
    if not np.all(np.diff(xp) > 0):
        raise ValueError("FE periodic interpolation grid is not strictly increasing.")

    # Reuse the identical exact-ZOH modal transfer for healthy and fault.
    # Retaining the signed force difference makes all mode and PR-source
    # contributions strictly additive at every sensor and time sample.
    n, ns, nm = time.size, sensor_local.size, frequency.size
    component = np.empty((n, ns, 3, nm), dtype=np.float64)
    healthy_response = np.zeros((n, ns), dtype=np.float64)
    delta_force = force_f - force_h
    dt = 1.0 / fs
    for im, hz in enumerate(frequency):
        fp = np.r_[shape_in[-1, im], shape_in[:, im], shape_in[0, im]]
        b, a = exact_zoh_acceleration_filter(2.0 * np.pi * hz, ZETA, dt)
        zi = lfilter_zi(b, a)
        for ip in range(3):
            psi = np.interp(theta_fe[:, ip], xp, fp)
            q_h = force_h[:, ip] * psi
            q_delta = delta_force[:, ip] * psi
            a_h, _ = lfilter(b, a, q_h, zi=zi * q_h[0])
            a_delta, _ = lfilter(b, a, q_delta, zi=zi * q_delta[0])
            sensor_mode = shape_sensor[:, im]
            healthy_response += a_h[:, None] * sensor_mode[None, :]
            component[:, :, ip, im] = a_delta[:, None] * sensor_mode[None, :]
    delta_response = component.sum(axis=(2, 3))
    fault_response = healthy_response + delta_response
    settled_start = int(round(SETTLED_START_S * fs))
    pr1_settled = delta_force[settled_start:-1, 0]
    fft_power = np.abs(np.fft.rfft((pr1_settled - np.mean(pr1_settled)) * np.hanning(pr1_settled.size))) ** 2
    fft_frequency = np.fft.rfftfreq(pr1_settled.size, d=dt)
    band_energy_fraction = {
        f"{low:g}-{high:g}_Hz": float(np.sum(fft_power[(fft_frequency >= low) & (fft_frequency < high)]) / np.sum(fft_power))
        for low, high in ((0.0, 5000.0), (5000.0, 15000.0), (15000.0, 20000.0), (20000.0, fs / 2.0 + 1.0))
    }

    # Independent whole-force replay checks decomposition, not just algebra.
    replay = np.zeros_like(delta_response)
    for im, hz in enumerate(frequency):
        fp = np.r_[shape_in[-1, im], shape_in[:, im], shape_in[0, im]]
        q_total = np.zeros(n)
        for ip in range(3):
            q_total += delta_force[:, ip] * np.interp(theta_fe[:, ip], xp, fp)
        b, a = exact_zoh_acceleration_filter(2.0 * np.pi * hz, ZETA, dt)
        a_total, _ = lfilter(b, a, q_total, zi=lfilter_zi(b, a) * q_total[0])
        replay += a_total[:, None] * shape_sensor[:, im][None, :]
    relative_reconstruction_error = np.linalg.norm(replay - delta_response) / np.linalg.norm(replay)
    if relative_reconstruction_error > 1e-8:
        raise ValueError(f"PR/mode decomposition does not reconstruct the FE response: {relative_reconstruction_error:g}")
    direct_asset = dict(asset)
    direct_asset["frequency_hz"] = frequency
    direct_asset["input_shape_normal"] = shape_in
    direct_asset["sensor_shape_radial"] = shape_sensor
    direct_response, _ = modal_zoh(time, delta_force, carrier, direct_asset, "corrected")
    direct_zoh_error = np.linalg.norm(direct_response - delta_response) / np.linalg.norm(direct_response)
    if direct_zoh_error > 1e-8:
        raise ValueError(f"Independent state-space exact-ZOH comparison failed: {direct_zoh_error:g}")

    events = event_segments(time, sp_loss, "SP") + event_segments(time, pr_loss, "PR")
    events.sort(key=lambda event: time[int(event["peak"])])
    if sum(event["kind"] == "SP" for event in events) != 5 or sum(event["kind"] == "PR" for event in events) != 5:
        raise ValueError("Reviewed one-second source-event counts changed.")

    metric_rows: list[dict[str, object]] = []
    detail_rows: list[dict[str, object]] = []
    event_rows: list[dict[str, object]] = []
    sensitivity_rows: list[dict[str, object]] = []
    event_winners: list[tuple[str, int, int]] = []
    for number, event in enumerate(events, start=1):
        first, last, peak = (int(event[key]) for key in ("first", "last", "peak"))
        window = event["window"]
        assert isinstance(window, slice)
        contact_global = float(np.rad2deg(theta_fe[peak, 0]) % 360.0)
        contact_local = float((sensor_global[0] - contact_global) % 360.0)
        distance = np.abs((contact_local - sensor_local + 180.0) % 360.0 - 180.0)
        nearest = int(np.argmin(distance))
        sensor_peaks = np.max(np.abs(delta_response[window]), axis=0)
        strongest = int(np.argmax(sensor_peaks))
        pr1_force_peak = float(np.max(np.abs(delta_force[window, 0])))
        three_pr_force_peak = float(np.max(np.linalg.norm(delta_force[window], axis=1)))
        event_winners.append((str(event["kind"]), nearest, strongest))
        event_rows.append({
            "event_number": number,
            "source_mesh": event["kind"],
            "fault_planet": 1,
            "source_start_s": float(time[first]),
            "source_peak_s": float(time[peak]),
            "source_end_s": float(time[last]),
            "pr_contact_angle_local_deg_at_source_peak": contact_local,
            "nearest_sensor_local_deg": float(sensor_local[nearest]),
            "strongest_sensor_local_deg": float(sensor_local[strongest]),
            "nearest_equals_strongest": int(nearest == strongest),
            "source_max_sp_loss_N_per_m": float(np.max(sp_loss[first : last + 1, 0])),
            "source_max_pr_loss_N_per_m": float(np.max(pr_loss[first : last + 1, 0])),
            "source_pr1_force_delta_peak_abs_N": pr1_force_peak,
            "source_three_pr_force_delta_vector_peak_N": three_pr_force_peak,
        })
        for isensor in range(ns):
            values = delta_response[window, isensor]
            local_peak_index = int(np.argmax(np.abs(values)))
            global_peak_index = window.start + local_peak_index
            abs_peak = float(np.abs(values[local_peak_index]))
            rms = float(np.sqrt(np.mean(values * values)))
            crossing = np.flatnonzero(np.abs(values) >= 0.5 * abs_peak)
            first50_index = window.start + int(crossing[0])
            total_peak = float(delta_response[global_peak_index, isensor])
            source_at_peak = component[global_peak_index, isensor].sum(axis=1)
            modal_at_peak = component[global_peak_index, isensor].sum(axis=0)
            metric_rows.append({
                "event_number": number,
                "source_mesh": event["kind"],
                "sensor_local_deg": float(sensor_local[isensor]),
                "pr_contact_angle_local_deg": contact_local,
                "contact_sensor_angular_distance_deg": float(distance[isensor]),
                "fault_increment_peak_abs_ms2": abs_peak,
                "fault_increment_window_rms_ms2": rms,
                "peak_per_pr1_force_delta_ms2_per_N": abs_peak / pr1_force_peak,
                "fault_increment_signed_at_peak_ms2": total_peak,
                "peak_time_s": float(time[global_peak_index]),
                "peak_lag_from_source_peak_ms": float((time[global_peak_index] - time[peak]) * 1000.0),
                "first_50pct_of_own_peak_time_s": float(time[first50_index]),
                "first_50pct_lag_from_source_start_ms": float((time[first50_index] - time[first]) * 1000.0),
                "active_pr_signed_share_at_total_peak": float(source_at_peak[0] / total_peak),
                "signed_cancellation_factor_at_total_peak": float(np.sum(np.abs(modal_at_peak)) / abs(total_peak)),
                "dominant_mode_index_at_total_peak": int(np.argmax(np.abs(modal_at_peak)) + 1),
                "dominant_mode_frequency_hz": float(frequency[np.argmax(np.abs(modal_at_peak))]),
            })
            for ip in range(3):
                for im in range(nm):
                    detail_rows.append({
                        "event_number": number,
                        "source_mesh": event["kind"],
                        "sensor_local_deg": float(sensor_local[isensor]),
                        "pr_input": ip + 1,
                        "mode_index": im + 1,
                        "mode_frequency_hz": float(frequency[im]),
                        "signed_component_at_total_peak_ms2": float(component[global_peak_index, isensor, ip, im]),
                        "total_signed_at_peak_ms2": total_peak,
                    })

    args.output.mkdir(parents=True, exist_ok=True)
    save_csv(args.output / "planet_fault_events.csv", event_rows)
    save_csv(args.output / "planet_sensor_event_metrics.csv", metric_rows)
    save_csv(args.output / "planet_pr_mode_peak_contributions.csv", detail_rows)
    for mode_count in (5, 11, 15, 20, 23):
        subset = component[:, :, :, :mode_count].sum(axis=(2, 3))
        for event_number, event in enumerate(events, start=1):
            window = event["window"]
            assert isinstance(window, slice)
            peak = int(event["peak"])
            contact_local = float((sensor_global[0] - np.rad2deg(theta_fe[peak, 0])) % 360.0)
            distances = np.abs((contact_local - sensor_local + 180.0) % 360.0 - 180.0)
            peak_by_sensor = np.max(np.abs(subset[window]), axis=0)
            first50 = []
            for isensor in range(ns):
                series = np.abs(subset[window, isensor])
                crossing = np.flatnonzero(series >= 0.5 * peak_by_sensor[isensor])
                first50.append(float(time[window.start + int(crossing[0])]))
            nearest = int(np.argmin(distances))
            strongest = int(np.argmax(peak_by_sensor))
            near_0_90 = 0 if distances[0] < distances[1] else 1
            early_0_90 = 0 if first50[0] < first50[1] else (1 if first50[1] < first50[0] else -1)
            sensitivity_rows.append({
                "mode_count": mode_count,
                "highest_mode_frequency_hz": float(frequency[mode_count - 1]),
                "event_number": event_number,
                "source_mesh": event["kind"],
                "local_pr_contact_deg": contact_local,
                "nearest_sensor_deg": float(sensor_local[nearest]),
                "strongest_sensor_deg": float(sensor_local[strongest]),
                "nearest_equals_strongest": int(nearest == strongest),
                "sensor_0_peak_abs_ms2": float(peak_by_sensor[0]),
                "sensor_90_peak_abs_ms2": float(peak_by_sensor[1]),
                "sensor_120_peak_abs_ms2": float(peak_by_sensor[2]),
                "sensor_240_peak_abs_ms2": float(peak_by_sensor[3]),
                "sensor_0_first50_s": first50[0],
                "sensor_90_first50_s": first50[1],
                "nearer_of_0_90_first50_earlier": int(near_0_90 == early_0_90),
            })
    save_csv(args.output / "planet_mode_count_sensitivity.csv", sensitivity_rows)
    signal_path = args.output / "planet_four_sensor_time_signals.csv"
    with signal_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_s"] + [f"healthy_sensor_{int(x):03d}deg_ms2" for x in sensor_local] + [f"planet_q090_sensor_{int(x):03d}deg_ms2" for x in sensor_local] + [f"planet_q090_minus_healthy_sensor_{int(x):03d}deg_ms2" for x in sensor_local])
        writer.writerows(np.column_stack((time, healthy_response, fault_response, delta_response)))

    summary: dict[str, object] = {
        "status": "exploratory: 23 of 60 FE modes under 25.6-kHz source Nyquist; 1 s settled segment",
        "healthy_source": str(healthy_path),
        "planet_q090_source": str(fault_path),
        "FE_asset": str(args.asset),
        "DOF": 18,
        "sampling_hz": fs,
        "source_duration_s": float(time[-1]),
        "settled_interval_s": [SETTLED_START_S, float(time[-1])],
        "modes_total": int(all_frequency.size),
        "modes_retained": int(nm),
        "retained_frequency_range_hz": [float(frequency.min()), float(frequency.max())],
        "excluded_modes_hz": all_frequency[~keep].tolist(),
        "modal_damping_ratio": ZETA,
        "angle_registration": "theta_FE = SENSOR_0_FE_angle - (carrier_angle + planet_phase), FE angles counter-clockwise; local angles clockwise from SENSOR_0",
        "event_window": "from source stiffness-loss start through 4 ms after source stiffness-loss end; all amplitudes use fault minus matched healthy FE response",
        "response_type": "radial sensor acceleration, signed, m/s^2; no normalization or amplification",
        "PR_input_force_convention": "mesh_force_pr_signed_N, matching pg_apply_fe_ring_transfer_causal's direct use of signed meshForceRingPlanet. Absolute polarity/phase needs FE unit-load and load-normal verification.",
        "mode_and_PR_reconstruction_relative_l2_error": float(relative_reconstruction_error),
        "independent_state_space_exact_ZOH_relative_l2_error": float(direct_zoh_error),
        "event_count_SP": sum(event["kind"] == "SP" for event in events),
        "event_count_PR": sum(event["kind"] == "PR" for event in events),
        "nearest_equals_strongest_by_mesh": {
            kind: {
                "count": sum(item[0] == kind for item in event_winners),
                "matches": sum(item[0] == kind and item[1] == item[2] for item in event_winners),
            }
            for kind in ("SP", "PR")
        },
        "mode_count_sensitivity_nearest_equals_strongest": {
            str(mode_count): sum(row["nearest_equals_strongest"] for row in sensitivity_rows if row["mode_count"] == mode_count)
            for mode_count in (5, 11, 15, 20, 23)
        },
        "mode_count_sensitivity_nearer_0_90_first50_earlier": {
            str(mode_count): sum(row["nearer_of_0_90_first50_earlier"] for row in sensitivity_rows if row["mode_count"] == mode_count)
            for mode_count in (5, 11, 15, 20, 23)
        },
        "source_PR1_force_increment_Hann_power_band_fraction": band_energy_fraction,
        "caveat": "The ANSYS modal mass/SI normalization has not been independently verified. This 23-mode, 51.2-kHz, one-second planet study is not quantitatively comparable to the earlier 60-mode, 131.072-kHz sun-crack study. First-50%-crossing and peak order are waveform metrics, not a measured travelling-wave speed.",
    }
    (args.output / "planet_path_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: summary[key] for key in ("modes_retained", "event_count_SP", "event_count_PR", "mode_and_PR_reconstruction_relative_l2_error", "nearest_equals_strongest_by_mesh")}, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--asset", type=Path, default=ASSET)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
