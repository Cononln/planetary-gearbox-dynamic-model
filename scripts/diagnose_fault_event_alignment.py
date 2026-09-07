"""Post-hoc envelope/signed-waveform check; reuse frozen shifts and event masks.

No fitting, candidate reassignment, parameter selection or SVD is done here.
The continuous-record Hilbert envelope avoids artificial short-window edges.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'runtime' / 'phase_python'))
sys.path.insert(0, str(ROOT / 'src_py'))
import numpy as np
from scipy import signal, ndimage
from event_alignment import cross_block_metrics, extract_events, interval


def save_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('result', type=Path)
    root = parser.parse_args().result.resolve()
    cfg = json.loads((root / 'runtime_config.json').read_text(encoding='utf-8'))
    selection = json.loads((root / 'selection.json').read_text(encoding='utf-8'))
    assert json.loads((root / 'evaluation_status.json').read_text(encoding='utf-8'))['complete']
    fs = cfg['working_rates_hz']['impulse']
    smooth_samples = 2 * int(round(.0005 * fs) // 2) + 1
    summary, strata = [], []
    original = {}
    with (root / 'metrics.csv').open(encoding='utf-8-sig') as stream:
        for row in csv.DictReader(stream):
            original[(row['condition'], row['method'])] = float(row['common_mask_phase'])
    for condition in ['BL', 'PF50']:
        folder = root / condition
        with np.load(folder / 'evaluated_events.npz') as archive:
            keys = ['t', 'psi', 'branch', 'common_mask', 'local_u_s']
            keys += [k for method in ['none', 'local_xcorr', 'constrained']
                     for k in [method, method + '_shift']]
            data = {key: archive[key] for key in keys}
        events = {key: data[key] for key in ['t', 'psi', 'branch']}
        wave = np.load(folder / f'band_{selection["band_index"]}' / 'waveform.npy', mmap_mode='r')
        envelope = np.empty(wave.shape, np.float32)
        for channel in range(3):
            analytic = signal.hilbert(wave[:, channel])
            envelope[:, channel] = ndimage.uniform_filter1d(np.abs(analytic), smooth_samples, mode='nearest')
        del analytic
        test_mask = np.zeros(len(data['t']), bool)
        for span in cfg['split_seconds']['test']:
            test_mask |= interval(data['t'], span)
        count = int((test_mask & data['common_mask']).sum())
        for method in ['none', 'local_xcorr', 'constrained']:
            env_windows, u = extract_events(envelope, fs, events, data[method + '_shift'])
            assert np.array_equal(u, data['local_u_s'])
            supports = []
            for representation, windows in [('signed', data[method]), ('envelope', env_windows)]:
                rows = cross_block_metrics(windows, events, data['common_mask'], cfg['split_seconds']['test'])
                assert rows and all(np.isfinite(row['phase']) for row in rows)
                supports.append({tuple(row[key] for key in ['branch', 'angle_bin', 'channel', 'block1', 'block2'])
                                 for row in rows})
                score = float(np.mean([row['phase'] for row in rows]))
                if representation == 'signed':
                    assert abs(score - original[(condition, method)]) < 1e-10
                summary.append(dict(condition=condition, method=method, representation=representation,
                                    coherence=score, matched_stratum_pairs=len(rows), test_events=count))
                strata.extend(dict(condition=condition, method=method, representation=representation, **row)
                               for row in rows)
            assert supports[0] == supports[1], 'Signed and envelope support mismatch'
            del env_windows
        del envelope, data
    save_csv(root / 'envelope_signed_diagnostic.csv', summary)
    save_csv(root / 'envelope_signed_strata.csv', strata)
    status = dict(scope='post_hoc_exploratory_diagnostic_of_completed_run',
                  refitting=False, new_event_selection=False, source_shifts='frozen evaluated_events arrays',
                  envelope='abs(Hilbert(continuous band waveform)), then centered moving average, then saved time shifts',
                  envelope_smoothing_samples=smooth_samples, envelope_smoothing_ms=smooth_samples/fs*1000,
                  scoring='same demeaned tapered unit-norm cross-block score as signed windows; before SVD',
                  mask='identical target event mask for all three methods and both representations',
                  signed_score_reproduced=True, representation_support_identical=True,
                  limits=['Envelope score is not a signed-phase metric or a fault-specificity proof.',
                          'No new wrong-order envelope controls or independent recordings were added.',
                          'Shifts can bring different source samples into the fixed core; no amplitude compensation is implied.'],
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  summaries=summary)
    (root / 'envelope_signed_diagnostic.json').write_text(json.dumps(status, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
