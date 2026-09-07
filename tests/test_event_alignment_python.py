"""Controlled shift recovery and leakage checks, distinct from measured evidence."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'runtime'/'phase_python'))
sys.path.insert(0, str(ROOT/'src_py'))
import json
import numpy as np
from event_alignment import path_design, fit_delay_model, predict_delays, sample_sinc, reconstruct_block


def run():
    rng = np.random.default_rng(4401)
    t = np.arange(.5, 16, 31/168/2)
    n = len(t)
    events = dict(t=t, psi=4*np.pi*t, branch=np.arange(n)%2)
    x = path_design(events['psi'], events['branch'])
    assert np.max(np.abs(x.sum(axis=1))) < 1e-12
    beta = np.zeros(20)
    beta[[2, 3, 4, 5, 12, 13, 14, 15]] = [.8, -.2, .3, .15, .65, .2, -.15, .2]
    tau = x@beta/1000
    common = .0015*np.sin(2*np.pi*t/3.5)
    truth = tau+common[:, None]
    offsets = np.stack([truth+rng.normal(0, .00012, truth.shape),
                        rng.uniform(-.008, .008, truth.shape),
                        rng.uniform(-.008, .008, truth.shape)], axis=2)
    scores = np.broadcast_to([7., 4., 3.5], offsets.shape).copy()
    missing = rng.random(truth.shape) < .06
    offsets[missing] = np.nan
    scores[missing] = -np.inf
    cand = dict(offsets=offsets, scores=scores)
    model = fit_delay_model(events, cand, [0, 8], smooth=.1)
    assert model['success'], model
    original = model['coeff'].copy()
    pred = predict_delays(events, cand, model, [[8, 16]])
    use = (t >= 8) & pred['accepted']
    before = np.sqrt(np.mean(truth[use]**2))
    after = np.sqrt(np.mean((pred['shift'][use]-truth[use])**2))
    assert after < .5*before and after < .0005, (before, after)
    a = predict_delays(events, cand, model, [[8, 16]], excluded=2)
    changed = {k: v.copy() for k, v in cand.items()}
    changed['offsets'][t >= 8, 2] = .007
    changed['scores'][t >= 8, 2] = 100
    b = predict_delays(events, changed, model, [[8, 16]], excluded=2)
    assert np.array_equal(a['shift'], b['shift'])
    assert np.array_equal(model['coeff'], original)
    fs = 51200
    grid = np.arange(fs)/fs
    y = np.cos(2*np.pi*6000*grid)+.4*np.sin(2*np.pi*10000*grid)
    q = np.arange(100, len(y)-100)+.37
    recovered = sample_sinc(y, q)
    expected = np.cos(2*np.pi*6000*q/fs)+.4*np.sin(2*np.pi*10000*q/fs)
    relative = np.linalg.norm(recovered-expected)/np.linalg.norm(expected)
    assert relative < .01, relative
    empty = dict(offsets=np.full_like(offsets, np.nan), scores=np.full_like(scores, -np.inf))
    rejected = fit_delay_model(events, empty, [0, 8])
    assert not rejected['success']
    ev = dict(t=np.array([.2, .4, .6, .8]))
    sig = np.tile(y[:, None], (1, 3))
    same, derivative, rejections = reconstruct_block(sig, fs, ev, np.zeros((4, 3)), np.ones(4, bool), [.1, .9])
    assert np.max(np.abs(same-sig[int(.1*fs):int(.9*fs)])) < 1e-6
    return dict(scope='Known-event candidate timing benchmark; not full vibration-event discovery or gearbox dynamics',
                passed=True, raw_delay_rmse_ms=before*1000, aligned_delay_rmse_ms=after*1000,
                test_event_coverage=float(use.sum()/np.sum(t >= 8)), interpolation_relative_error=relative,
                heldout_channel_invariance=True, path_coefficients_frozen=True, zero_evidence_rejected=True,
                zero_shift_reconstruction_max_error=float(np.max(np.abs(same-sig[int(.1*fs):int(.9*fs)]))))


if __name__ == '__main__':
    result = run()
    path = ROOT/'results'/'event_alignment_python_tests_20260906.json'
    path.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2), flush=True)
