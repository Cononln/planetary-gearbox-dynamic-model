"""Recheck persisted signed-waveform outputs without using the metric helper."""
from pathlib import Path
import csv
import json
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'runtime'/'phase_python'))
import numpy as np


def main():
    out = Path(sys.argv[1]).resolve()
    with (out/'metrics.csv').open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle))
    max_fft_error, max_metric_error, max_fusion_error = 0., 0., 0.
    for condition in ['BL', 'PF50']:
        with np.load(out/f'{condition}_signed_results.npz') as data:
            f = data['frequency_hz']
            bins = np.arange(len(f))
            distance = abs(bins-np.rint(bins/2604)*2604)
            centers = bins[(f >= 4000) & (f <= 10000) & (bins % 84 == 0) & (distance > 103)]
            target = centers[:, None]+np.arange(-1, 2)
            background = centers[:, None]+np.r_[np.arange(-10, -4), np.arange(5, 11)]
            for method in ['none', 'whole_waveform', 'tidal_waveform']:
                for rank in [0, 1, 2]:
                    prefix = f'{method}_rank{rank}'
                    x = data[prefix+'_channels'].astype(float)
                    fused = (x-x.mean(axis=1, keepdims=True)).mean(axis=0)
                    fusion_error = np.max(abs(fused-data[prefix+'_fused']))/max(np.max(abs(fused)), 1e-30)
                    max_fusion_error = max(max_fusion_error, float(fusion_error))
                    z = np.fft.rfft(fused)
                    amp = 2*abs(z)/len(fused)
                    amp[0] /= 2
                    amp[-1] /= 2
                    fft_error = np.max(abs(amp-data[prefix+'_amplitude']))/max(np.max(amp), 1e-30)
                    max_fft_error = max(max_fft_error, float(fft_error))
                    computed = 10*np.log10(np.mean(amp[target]**2)/np.mean(amp[background]**2))
                    saved = next(r for r in rows if r['condition'] == condition and r['method'] == method and int(r['rank']) == rank)
                    max_metric_error = max(max_metric_error, abs(float(saved['fault_spaced_comb_background_db'])-computed))
                    assert np.isclose(np.sqrt(np.mean(fused*fused)), float(saved['signed_fused_rms']), rtol=1e-12)
    assert max_fft_error < 1e-6 and max_fusion_error < 1e-6 and max_metric_error < 1e-10
    audit = dict(passed=True, evaluated_outputs=18, max_relative_saved_fft_error=max_fft_error,
                 max_relative_saved_fusion_error=max_fusion_error,
                 max_recomputed_comb_difference_db=float(max_metric_error))
    (out/'persisted_output_audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
