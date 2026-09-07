"""Known-arrival waveform controls; not a gearbox dynamics or fault classifier test."""
from pathlib import Path
import json
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'runtime'/'phase_python'))
sys.path.insert(0, str(ROOT/'src_py'))
import numpy as np
from event_arrival import (envelope_response, onset_markers, arrival_templates, envelope_candidates,
                           local_result, fit_relative_path, joint_result)
from event_alignment import path_design, interval


def benchmark(seed):
    rng = np.random.default_rng(seed); fs = 51200
    t = np.arange(.3, 11.7, 31/168/2)
    events = dict(t=t, psi=4*np.pi*t, branch=np.arange(len(t))%2)
    beta = np.zeros(20); beta[[2, 3, 4, 5, 12, 13, 14, 15]] = [.7, .2, .1, -.1, .6, .1, .15, -.1]
    truth = (path_design(events['psi'], events['branch'])@beta)/1000
    truth += (.002*np.sin(.8*t)+.0006*np.sin(3*t))[:, None]
    grid = np.arange(fs*12)/fs
    wave = rng.normal(0, .025, (len(grid), 3))
    for n, t0 in enumerate(t):
        for j in range(3):
            center = t0+truth[n, j]
            ids = np.flatnonzero((grid >= center) & (grid < center+.012))
            u = grid[ids]-center
            pulse = (1-np.exp(-u/.00018))*np.exp(-u/.0013)
            wave[ids, j] += rng.uniform(.8, 1.2)*pulse*np.cos(2*np.pi*[5300, 6700, 8200][j]*u+rng.uniform(0, 2*np.pi))
    env = envelope_response(wave.astype(np.float32), fs)
    marker = onset_markers(env, fs, events)
    templates, _ = arrival_templates(env, fs, events, marker, [0, 5])
    cand = envelope_candidates(env, fs, events, templates)
    local = local_result(cand); model = fit_relative_path(events, cand, [0, 5])
    joint = joint_result(events, cand, model, [[5, 12]], .5)
    use = interval(t, [5, 12]) & local['accepted'] & joint['accepted']
    rows = []
    for method, shift in [('none', np.zeros_like(truth)), ('envelope_xcorr', local['shift']), ('joint_envelope', joint['shift'])]:
        error = (shift[use]-truth[use])*1000
        rows.append(dict(seed=seed, method=method, timing_mae_ms=float(np.mean(np.abs(error))),
                         timing_rmse_ms=float(np.sqrt(np.mean(error**2))), covered_events=int(use.sum()),
                         scheduled_events=int(interval(t, [5, 12]).sum())))
    assert use.sum() > .9*interval(t, [5, 12]).sum()
    assert rows[1]['timing_mae_ms'] < .4 and rows[2]['timing_mae_ms'] < .5, rows
    original = model['coeff'].copy()
    pred = joint_result(events, cand, model, [[5, 12]], .5, excluded=2)
    altered = {k: value.copy() for k, value in cand.items()}
    altered['offsets'][t >= 5, 2] = .007
    altered['scores'][t >= 5, 2] = .999
    pred2 = joint_result(events, altered, model, [[5, 12]], .5, excluded=2)
    assert np.array_equal(pred['predicted'], pred2['predicted'])
    assert np.array_equal(pred['accepted'], pred2['accepted'])
    assert np.array_equal(original, model['coeff'])
    empty = dict(offsets=np.full_like(cand['offsets'], np.nan), scores=np.full_like(cand['scores'], -np.inf))
    reject = joint_result(events, empty, model, [[5, 12]], .5)
    assert not reject['accepted'].any() and not reject['shift'].any()
    assert np.max(np.abs(pred['path'].sum(axis=1))) < 1e-12
    return rows


if __name__ == '__main__':
    rows = sum([benchmark(seed) for seed in [4401, 4402, 4403]], [])
    result = dict(passed=True, scope='Known scheduled arrivals in synthetic band waveforms with random ringing phase; not full unsupervised fault discovery',
                  heldout_channel_invariance=True, path_frozen=True, zero_evidence_rejected=True,
                  common_reference_zero_mean=True, results=rows)
    (ROOT/'results'/'event_arrival_v2_tests.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2), flush=True)
