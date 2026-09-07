"""Independent arithmetic checks of persisted sun-fault trial outputs."""
from pathlib import Path
import csv
import json
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'runtime'/'phase_python'))
import numpy as np


def main():
    out = Path(sys.argv[1]).resolve()
    manifest = json.loads((out/'run_manifest.json').read_text(encoding='utf-8'))
    with (out/'metrics.csv').open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle))
    with (out/'periods.csv').open(encoding='utf-8-sig', newline='') as handle:
        periods = list(csv.DictReader(handle))
    assert len(periods) == manifest['rows'] == 111
    assert manifest['kinematics'] == dict(carrier_turns=1, mesh_cycles=84, fault_cycles=12)
    for left, right in zip(periods[:-1], periods[1:]):
        assert float(left['local_end_s']) == float(right['local_start_s'])
    relative_error, db_error = [], []
    with np.load(out/'Sun050_signed_results.npz') as data:
        f = data['frequency_hz']
        bins = np.arange(len(f))
        centers = bins[(f >= 4000) & (f <= 10000) & (bins % 12 == 0)
                       & (abs(bins-np.rint(bins/84)*84) > 13)]
        target = centers[:, None]+np.arange(-1, 2)
        background = centers[:, None]+np.r_[np.arange(-10, -4), np.arange(5, 11)]
        for method in ['none', 'whole_waveform', 'tidal_waveform']:
            lag = data[method+'_lags_ms']
            assert lag.shape == (111, 3) and np.all(lag[:, 0] == 0)
            if method != 'tidal_waveform':
                assert np.all(lag == lag[0])
            for rank in [0, 1, 2]:
                key = f'{method}_rank{rank}'
                x = data[key+'_channels'].astype(float)
                assert x.shape == (3, 16384)
                centered = x-x.mean(axis=1, keepdims=True)
                fused = centered.mean(axis=0)
                amp = 2*abs(np.fft.rfft(fused))/len(fused)
                amp[[0, -1]] /= 2
                relative_error.append(float(np.max(abs(amp-data[key+'_amplitude']))/np.max(amp)))
                relative_error.append(float(np.max(abs(fused-data[key+'_fused']))/np.max(abs(fused))))
                ratio = 10*np.log10(np.mean(amp[target]**2)/np.mean(amp[background]**2))
                saved = next(r for r in rows if r['method'] == method and int(r['rank']) == rank)
                db_error.append(abs(ratio-float(saved['fault_spaced_comb_background_db'])))
                assert np.isclose(np.sqrt(np.mean(fused*fused)), float(saved['signed_fused_rms']), rtol=1e-12)
                if method != 'none':
                    assert np.array_equal(x[0], data[f'none_rank{rank}_channels'][0])
    assert max(relative_error) < 1e-6 and max(db_error) < 1e-10
    audit = dict(passed=True, outputs_checked=9, full_tidal_rows_verified=111,
                 max_saved_relative_error=max(relative_error), max_comb_recompute_difference_db=float(max(db_error)),
                 same_spans=True, s1_fixed=True, whole_lags_constant=True)
    (out/'persisted_output_audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
