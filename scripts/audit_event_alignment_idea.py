"""Deterministic algebra checks for the event-alignment research proposal.

These are counterexamples and identities, not gearbox simulation or real-data
validation. Uses only the Python standard library. Run from any directory:
    py -3.14 scripts/audit_event_alignment_idea.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path


def dot(a, b):
    return math.fsum(x * y for x, y in zip(a, b, strict=True))


def mean_rows(rows):
    return [math.fsum(column) / len(rows) for column in zip(*rows, strict=True)]


def centered(row):
    average = math.fsum(row) / len(row)
    return [x - average for x in row]


def unit(row):
    values = centered(row)
    norm = math.sqrt(dot(values, values))
    if norm <= 0:
        raise ValueError("A zero-energy event has no normalized phase score.")
    return [x / norm for x in values]


def old_c_tsa(rows):
    q = [unit(row) for row in rows]
    average = mean_rows(q)
    count = len(q)
    return (count * dot(average, average) - 1) / (count - 1)


def mean_pair_correlation(rows):
    q = [unit(row) for row in rows]
    values = [dot(q[i], q[j]) for i in range(len(q)) for j in range(i)]
    return math.fsum(values) / len(values)


def cross_block_fraction(blocks):
    """Signed cross-block coherent energy / mean individual event energy.

    Events must already use the same fixed waveform support, channel, event
    type and compatible carrier-angle stratum. This function performs NO
    alignment, event selection, sign flipping, gain fitting or normalization.
    """
    if len(blocks) < 2 or any(not block for block in blocks):
        raise ValueError("At least two nonempty held-out blocks are required.")
    width = len(blocks[0][0])
    if width == 0 or any(len(row) != width for block in blocks for row in block):
        raise ValueError("Every event must have the same nonempty support.")
    if any(not math.isfinite(x) for block in blocks for row in block for x in row):
        raise ValueError("Events must contain finite samples.")
    averages = [mean_rows(block) for block in blocks]
    count = len(blocks)
    cross_energy = 2 * math.fsum(
        dot(averages[i], averages[j])
        for i in range(count) for j in range(i)
    ) / (count * (count - 1))
    individual_energy = math.fsum(
        math.fsum(dot(row, row) for row in block) / len(block)
        for block in blocks
    ) / count
    if individual_energy <= 0:
        raise ValueError("Zero energy is undefined, not perfect alignment.")
    return cross_energy / individual_energy


def scaled(rows, factor):
    return [[factor * x for x in row] for row in rows]


def phase_cross_block_score(blocks):
    """Normalize the scoring copy only; preserve the sign of each event."""
    return cross_block_fraction([[unit(row) for row in block] for block in blocks])


def run_audit():
    checks = []

    def record(name, values, passed, meaning):
        if not passed:
            raise AssertionError(f"{name}: {values}")
        checks.append(dict(name=name, values=values, passed=True, meaning=meaning))

    rows = [[math.sin((n + 1) * 0.13 + j * 0.2) +
             0.2 * math.cos((n + 1) * (j + 1) * 0.17)
             for n in range(80)] for j in range(9)]
    c_tsa, mean_corr = old_c_tsa(rows), mean_pair_correlation(rows)
    record("C_TSA_is_mean_pair_correlation",
           {"C_TSA": c_tsa, "mean_pair_correlation": mean_corr},
           abs(c_tsa - mean_corr) < 1e-12,
           "The proposed old C_TSA is exactly average Pearson correlation "
           "after per-event demeaning and unit-norm scaling.")

    mu = mean_rows(rows)
    total = math.fsum(dot(row, row) for row in rows) / len(rows)
    coherent = dot(mu, mu)
    variance = total - coherent  # population variance, divisor N
    pcc_squared = coherent / variance
    ratio = coherent / total
    record("same_sample_energy_ratio_is_PCC_transform",
           {"energy_ratio": ratio, "PCC2_over_1_plus_PCC2":
            pcc_squared / (1 + pcc_squared)},
           abs(ratio - pcc_squared / (1 + pcc_squared)) < 1e-12,
           "Under population-STD convention, ordinary TSA energy fraction "
           "is PCC squared / (1 + PCC squared); it is not a new statistic.")

    pulse = unit([0.0, 0.0, 1.0, -0.7, 0.3, -0.1, 0.0, 0.0])
    negative = [-x for x in pulse]
    opposing = [pulse] * 4 + [negative] * 4
    average = mean_rows(opposing)
    record("rank_one_does_not_imply_TSA_retention",
           {"rank": 1, "TSA_energy": dot(average, average),
            "C_TSA": old_c_tsa(opposing)},
           dot(average, average) < 1e-14 and old_c_tsa(opposing) < 0,
           "Rows +s and -s have rank one, but their TSA cancels exactly.")

    opposed_blocks = [[pulse, pulse], [negative, negative]]
    restored_blocks = [[pulse, pulse], [pulse, pulse]]
    r_opposed = cross_block_fraction(opposed_blocks)
    r_restored = cross_block_fraction(restored_blocks)
    record("center_alignment_misses_polarity",
           {"opposed_phase_R": r_opposed, "same_phase_R": r_restored},
           abs(r_opposed + 1) < 1e-12 and abs(r_restored - 1) < 1e-12,
           "Identical envelope arrival times can coexist with opposite "
           "signed ringing phases. A peak-timing metric cannot detect this.")

    background = [[pulse, pulse], [pulse, pulse]]
    r_background = cross_block_fraction(background)
    specificity = r_background - cross_block_fraction(background)
    record("repeated_background_requires_negative_controls",
           {"raw_R": r_background, "matched_control_corrected_S": specificity},
           abs(r_background - 1) < 1e-12 and abs(specificity) < 1e-12,
           "A repeated healthy waveform also scores one. Matched controls "
           "are needed; phase repeatability alone is not fault specificity.")

    # Orthogonal noise makes a transparent deterministic energy calculation.
    v = unit([0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 1.0, -1.0])
    projection = dot(v, pulse)
    noise = unit([a - projection * b for a, b in zip(v, pulse, strict=True)])
    positive_noise = [s + 0.5 * e for s, e in zip(pulse, noise, strict=True)]
    negative_noise = [s - 0.5 * e for s, e in zip(pulse, noise, strict=True)]
    fault_blocks = [[positive_noise, negative_noise]] * 2
    control_blocks = [[noise, [-x for x in noise]]] * 2
    r_fault = cross_block_fraction(fault_blocks)
    r_control = cross_block_fraction(control_blocks)
    tenfold = [scaled(block, 10.0) for block in fault_blocks]
    record("coherent_energy_and_gain_invariance",
           {"fault_R": r_fault, "control_R": r_control,
            "S": r_fault - r_control,
            "tenfold_sensitivity_R": cross_block_fraction(tenfold)},
           abs(r_fault - 0.8) < 1e-12 and abs(r_control) < 1e-12 and
           abs(r_fault - cross_block_fraction(tenfold)) < 1e-12,
           "R equals coherent energy fraction in this controlled case and "
           "is invariant to one fixed channel gain; physical amplitude "
           "comparisons across unknown gains remain unsupported.")

    amplitude_varied = [[pulse, [2 * x for x in pulse]],
                        [[3 * x for x in pulse], [4 * x for x in pulse]]]
    r_energy = cross_block_fraction(amplitude_varied)
    r_phase = phase_cross_block_score(amplitude_varied)
    record("phase_score_separates_amplitude_variation",
           {"raw_energy_R": r_energy, "normalized_phase_R": r_phase},
           abs(r_energy - 0.7) < 1e-12 and abs(r_phase - 1) < 1e-12,
           "Identical signed waveforms with different positive amplitudes "
           "are phase aligned. Normalize only the scoring copy to remove "
           "this amplitude dependence, and report raw-energy R separately.")

    shuffled = [list(reversed(block)) for block in fault_blocks]
    record("shuffling_event_rows_is_not_a_phase_null",
           {"original_R": r_fault, "row_shuffled_R": cross_block_fraction(shuffled)},
           abs(r_fault - cross_block_fraction(shuffled)) < 1e-12,
           "Changing row order preserves TSA and is not a valid "
           "wrong-phase negative control.")

    psi = [2 * math.pi * i / 97 for i in range(97)]
    alpha = [0.0, 2 * math.pi / 3, 4 * math.pi / 3]
    differential = {}
    for order in (1, 3, 6):
        maximum = 0.0
        for angle in psi:
            h = [math.cos(order * (angle - a)) for a in alpha]
            h_mean = math.fsum(h) / 3
            maximum = max(maximum, max(abs(x - h_mean) for x in h))
        differential[str(order)] = maximum
    record("three_sensor_common_3k_blind_modes",
           differential,
           differential["1"] > 0.5 and differential["3"] < 1e-12 and
           differential["6"] < 1e-12,
           "For a shifted common path function, 3k angular orders are "
           "identical at sensors 120 degrees apart and vanish on "
           "channel-mean subtraction. They are not identifiable relative delays.")

    max_gauge_error = 0.0
    for angle in psi:
        drift = 0.002 * math.sin(0.4 * angle)
        gauge = 0.001 * math.cos(3 * angle)
        for sensor_angle in alpha:
            delay = 0.0005 * math.cos(angle - sensor_angle)
            original = drift + delay
            transformed = (drift + gauge) + (delay - gauge)
            max_gauge_error = max(max_gauge_error, abs(original - transformed))
    record("common_drift_path_gauge",
           {"max_observation_difference_seconds": max_gauge_error},
           max_gauge_error < 1e-14,
           "Adding an arbitrary common function to drift and subtracting "
           "it from all paths changes no arrival observation. Gauge fixing "
           "does not prove physical speed-path separation.")

    zs, zp, zr = 21, 31, 84
    fsun = 600 / 60
    fc = fsun * zs / (zs + zr)
    fm = zr * fc
    repeat_turns = zp // math.gcd(zp, zr)
    repeat_sec = repeat_turns / fc
    record("nominal_event_recurrence",
           {"carrier_hz": fc, "mesh_hz": fm,
            "planet_single_branch_hz": fm / zp,
            "planet_two_branch_rate_hz": 2 * fm / zp,
            "sun_all_three_planets_event_rate_hz": 3 * fm / zs,
            "planet_same_carrier_position_repeat_turns": repeat_turns,
            "repeat_seconds": repeat_sec,
            "complete_repeats_in_60s_without_guards": math.floor(60 / repeat_sec)},
           fm == 168 and repeat_turns == 31 and repeat_sec == 15.5,
           "Conditional on current code tooth counts, exact planet "
           "event/carrier-position recurrence is 31 carrier turns. Do not "
           "assume every tidal period has twelve planet fault events.")

    try:
        cross_block_fraction([[[0.0, 0.0]], [[0.0, 0.0]]])
    except ValueError:
        zero_rejected = True
    else:
        zero_rejected = False
    record("undefined_zero_energy_is_rejected", {}, zero_rejected,
           "No signal is not evidence of perfect alignment.")

    return {
        "scope": "Deterministic mathematical audit; no gearbox or measured-data experiment.",
        "all_passed": all(row["passed"] for row in checks),
        "check_count": len(checks),
        "checks": checks,
        "does_not_establish": [
            "New alignment algorithm performance on measured records",
            "Fault specificity under arbitrary unmatched background controls",
            "Independence of held-out blocks under adaptive same-window shift fitting",
            "Absolute structural transfer-path phase or speed-path identifiability",
            "Literature-wide novelty or journal acceptance",
        ],
    }


if __name__ == "__main__":
    result = run_audit()
    project = Path(__file__).resolve().parents[1]
    output = project / "results" / "event_alignment_idea_audit_20260906.json"
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
