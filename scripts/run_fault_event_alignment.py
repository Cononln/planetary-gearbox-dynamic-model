"""Run the actual BL/PF50 pilot with frozen-development settings and controls.

Use --run to acknowledge measured processing; --stage prepare caches only
development inference. --resume PATH evaluates a completed preparation cache.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
import shutil
from datetime import datetime
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'runtime'/'phase_python'))
sys.path.insert(0, str(ROOT/'src_py'))
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'runtime'/'matplotlib_feasibility'))
os.environ.setdefault('OMP_NUM_THREADS', '2')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
import numpy as np
from scipy import signal
from event_alignment import *


def jsonable(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    raise TypeError(type(obj).__name__)


def save_json(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, default=jsonable)+'\n', encoding='utf-8')


def save_csv(path, rows):
    if not rows:
        path.write_text('status\nno_valid_observations\n', encoding='utf-8')
        return
    with path.open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def spans(cfg):
    s = cfg['split_seconds']
    return [s['calibration'], s['validation']]+s['test']


def validation_score(events, cand, model, cfg):
    val = cfg['split_seconds']['validation']
    available = interval(events['t'], val)
    errors, count = [], 0
    for j in range(3):
        pred = predict_delays(events, cand, model, [val], excluded=j)
        # Strongest precomputed evidence candidate is independent of prediction.
        observed = cand['offsets'][:, j, 0]
        use = available & pred['accepted'] & np.isfinite(observed)
        errors.extend(np.abs(pred['shift'][use, j]-observed[use])*1000)
        count += int(use.sum())
    coverage = count/max(3*available.sum(), 1)
    mae = float(np.mean(errors)) if errors else float('inf')
    return mae+5*(1-coverage), mae, coverage


def prepare(out, cfg, reuse=None):
    rows = []
    for rec in cfg['records']:
        condition = rec['condition']
        folder = out/condition
        folder.mkdir()
        print(f'{condition}: loading source and deriving vibration motion', flush=True)
        if reuse is None:
            high, low, source_fs, encoder_edges = load_record(rec['path'], cfg)
            motion = estimate_motion(low, cfg)
            del low
            np.savez(folder/'motion.npz', **motion)
            np.savez(folder/'encoder_evaluation_only.npz', a=encoder_edges[0], b=encoder_edges[1])
            save_json(folder/'source.json', dict(source_fs=source_fs, samples_hi=len(high),
                      source_bytes=Path(rec['path']).stat().st_size, source_mtime_ns=Path(rec['path']).stat().st_mtime_ns,
                      vibration_columns_zero_based=cfg['vibration_columns']))
        else:
            for name in ['motion.npz', 'encoder_evaluation_only.npz', 'source.json']:
                shutil.copy2(reuse/condition/name, folder/name)
            motion = load_npz(folder/'motion.npz')
        print(f'{condition}: mesh {motion["fm"]:.6f} Hz; mean-frequency anchor {motion["slope_bias_hz"]:.6f} Hz', flush=True)
        for bi, band in enumerate(cfg['transient']['candidate_bands_hz']):
            bd = folder/f'band_{bi}'
            bd.mkdir()
            if reuse is None:
                waveform, initial_d, energy = make_evidence(high, motion, band, cfg)
                events = schedule(initial_d, motion, cfg)
                d = contrast_from_energy(energy, motion, cfg, events)
                cand = candidates(events, d, motion, cfg)
                np.save(bd/'waveform.npy', waveform)
                np.save(bd/'evidence.npy', d)
                np.save(bd/'energy.npy', energy)
                np.savez(bd/'events.npz', **events)
                np.savez(bd/'candidates.npz', **cand)
                del waveform, energy, initial_d, d
            else:
                for name in ['waveform.npy', 'evidence.npy', 'energy.npy', 'events.npz', 'candidates.npz']:
                    shutil.copy2(reuse/condition/f'band_{bi}'/name, bd/name)
                events = load_npz(bd/'events.npz'); cand = load_npz(bd/'candidates.npz')
            detected = np.isfinite(cand['offsets'][:, :, 0]).sum(axis=1) >= 2
            for si, smooth in enumerate(cfg['delay_model']['common_smoothness_candidates']):
                model = fit_delay_model(events, cand, cfg['split_seconds']['calibration'], smooth=smooth)
                loss, mae, coverage = validation_score(events, cand, model, cfg)
                save_json(bd/f'model_{si}.json', model)
                row = dict(condition=condition, band_index=bi, band_low_hz=band[0], band_high_hz=band[1],
                           smooth_index=si, smooth=smooth, model_success=model['success'],
                           detected_fraction=float(detected.mean()), validation_mae_ms=mae,
                           validation_coverage=coverage, selection_loss=loss)
                rows.append(row)
                print(f'{condition} band {band} lambda={smooth}: detected={detected.mean():.3f}, '
                      f'validation MAE={mae:.3f} ms coverage={coverage:.3f} success={model["success"]}', flush=True)
            save_csv(out/'development_candidates.csv', rows)
        if reuse is None:
            del high
    selections = {}
    totals = []
    for bi in range(len(cfg['transient']['candidate_bands_hz'])):
        selected = {}
        for rec in cfg['records']:
            condition = rec['condition']
            candidates_ = [r for r in rows if r['condition'] == condition and r['band_index'] == bi]
            selected[condition] = min(candidates_, key=lambda r: (r['selection_loss'], r['smooth_index']))
        totals.append((float(np.mean([r['selection_loss'] for r in selected.values()])), bi, selected))
    loss, bi, selections = min(totals, key=lambda r: (r[0], r[1]))
    if not np.isfinite(loss):
        bi = cfg['transient']['candidate_bands_hz'].index(cfg['transient']['initial_band_hz'])
        selections = {rec['condition']: min([r for r in rows if r['condition'] == rec['condition'] and r['band_index'] == bi],
                                          key=lambda r: r['smooth_index']) for rec in cfg['records']}
    save_json(out/'selection.json', dict(band_index=bi, common_band_hz=cfg['transient']['candidate_bands_hz'][bi],
                                       selections=selections, selection_data='calibration_and_validation_only',
                                       selected_models_available=bool(np.isfinite(loss))))
    print(f'Frozen common band: {cfg["transient"]["candidate_bands_hz"][bi]}', flush=True)


def load_npz(path):
    with np.load(path) as data:
        return {k: data[k] for k in data.files}


def global_xcorr(wave, fs, events, motion, cfg, accepted):
    # One scalar shift per relative-planet rotation, not a full tidal repeat.
    tp = cfg['gear_teeth']['planet']/float(motion['fm'])
    a_idx = np.flatnonzero(events['branch'] == 0)
    sub = {k: v[a_idx] for k, v in events.items() if isinstance(v, np.ndarray) and v.shape == events['t'].shape}
    z = np.zeros((len(a_idx), 3))
    cyc, _ = extract_events(wave, fs, sub, z, core=(-.15*tp, .85*tp))
    templates = make_templates(cyc, sub, cfg['split_seconds']['calibration'])
    binned = (np.mod(sub['psi'], 2*np.pi)/(np.pi/3)).astype(int)
    shifts = np.zeros((len(events['t']), 3))
    search = round(.008*fs)
    width = cyc.shape[2]
    for n, (t0, bin_) in enumerate(zip(sub['t'], binned)):
        center = round((t0-.15*tp)*fs)
        group = (events['t'] >= t0-.00001) & (events['t'] < t0+tp-.00001) & accepted
        for j in range(3):
            template = templates.get((0, int(bin_), j))
            if template is None or np.linalg.norm(template) < 1e-20:
                continue
            source = wave[center-search:center+width+search, j].astype(float)
            score = signal.correlate(source, template, mode='valid', method='fft')
            sums = np.r_[0, np.cumsum(source)]; squares = np.r_[0, np.cumsum(source**2)]
            variance = squares[width:]-squares[:-width]-(sums[width:]-sums[:-width])**2/width
            score /= np.maximum(np.sqrt(np.maximum(variance, 0))*np.linalg.norm(template), 1e-20)
            shifts[group, j] = (int(np.argmax(score))-search)/fs
    return shifts


def add_scores(windows, events, mask, cfg, label, condition, hypothesis, target):
    rows = cross_block_metrics(windows, events, mask, cfg['split_seconds']['test'])
    for row in rows:
        row.update(condition=condition, hypothesis=hypothesis, method=label)
    target.extend(rows)
    return float(np.mean([r['phase'] for r in rows])) if rows else float('nan')


def evaluate(out, cfg):
    selection = json.loads((out/'selection.json').read_text(encoding='utf-8'))
    bi = selection['band_index']
    metrics, phase_rows, spectrum_rows, arrival_rows, coverage_rows = [], [], [], [], []
    for rec in cfg['records']:
        condition = rec['condition']; folder = out/condition; bd = folder/f'band_{bi}'
        motion = load_npz(folder/'motion.npz')
        wave = np.load(bd/'waveform.npy', mmap_mode='r')
        energy = np.load(bd/'energy.npy', mmap_mode='r')
        events = load_npz(bd/'events.npz'); cand = load_npz(bd/'candidates.npz')
        si = selection['selections'][condition]['smooth_index']
        model = json.loads((bd/f'model_{si}.json').read_text(encoding='utf-8')); model['coeff'] = np.array(model['coeff'])
        fs = cfg['working_rates_hz']['impulse']
        allspans = spans(cfg)
        accepted = np.isfinite(cand['offsets'][:, :, 0]).sum(axis=1) >= 2
        raw, u = extract_events(wave, fs, events, np.zeros((len(events['t']), 3)))
        templates = make_templates(raw, events, cfg['split_seconds']['calibration'])
        print(f'{condition}: fitting operational comparisons and held-out-channel inference', flush=True)
        proposed = predict_delays(events, cand, model, allspans)
        local = local_xcorr(wave, fs, events, templates, accepted)
        whole = global_xcorr(wave, fs, events, motion, cfg, accepted)
        methods = dict(none=dict(shift=np.zeros_like(local), accepted=accepted),
                       global_xcorr=dict(shift=whole, accepted=accepted),
                       local_xcorr=dict(shift=local, accepted=accepted), constrained=proposed)
        for variant in ['independent', 'no_path_smoothing']:
            alt = fit_delay_model(events, cand, cfg['split_seconds']['calibration'], model['smooth'], variant)
            save_json(folder/f'{variant}_model.json', alt)
            methods[variant] = predict_delays(events, cand, alt, allspans)
        held = dict(shift=np.zeros_like(local), accepted=np.ones(len(local), bool))
        for j in range(3):
            pred = predict_delays(events, cand, model, allspans, excluded=j)
            held['shift'][:, j] = pred['shift'][:, j]
            held['accepted'] &= pred['accepted']
            observed = cand['offsets'][:, j, 0]
            for block, span in enumerate(cfg['split_seconds']['test']):
                use = interval(events['t'], span) & pred['accepted'] & np.isfinite(observed)
                arrival_rows.append(dict(condition=condition, channel=j+1, block=block,
                    event_count=int(use.sum()), coverage=float(use.sum()/max(interval(events['t'], span).sum(), 1)),
                    uncorrected_proxy_mae_ms=float(np.mean(np.abs(observed[use]))*1000) if use.any() else np.nan,
                    predicted_proxy_mae_ms=float(np.mean(np.abs(pred['shift'][use, j]-observed[use]))*1000) if use.any() else np.nan))
        methods['heldout_constrained'] = held
        common_mask = accepted & proposed['accepted']
        saved = dict(raw=raw, local_u_s=u, t=events['t'], psi=events['psi'], branch=events['branch'],
                     common_mask=common_mask, path_s=proposed['path'], common_shift_s=proposed['common'])
        for label, result in methods.items():
            print(f'{condition}: {label}, accepted {result["accepted"].mean():.3f}', flush=True)
            windows = raw if label == 'none' else extract_events(wave, fs, events, result['shift'])[0]
            scoring_mask = common_mask if label in ['none', 'global_xcorr', 'local_xcorr', 'constrained'] else common_mask & result['accepted']
            phase = add_scores(windows, events, scoring_mask, cfg, label, condition, 'target', phase_rows)
            denoised, svd_groups = frozen_svd(windows, events, result['accepted'], cfg)
            saved[label] = windows
            saved[label+'_svd'] = denoised
            saved[label+'_shift'] = result['shift']
            saved[label+'_mask'] = result['accepted']
            for block, span in enumerate(cfg['split_seconds']['test']):
                use = interval(events['t'], span)
                coverage = float(result['accepted'][use].mean())
                coverage_rows.append(dict(condition=condition, method=label, block=block, total_events=int(use.sum()),
                                          accepted_events=int(np.sum(use & result['accepted'])), coverage=coverage))
                yy, mind, reject = reconstruct_block(wave, fs, events, result['shift'], result['accepted'], span)
                order, amplitude, scores = family_spectrum(yy, fs, span, motion)
                metrics.append(dict(condition=condition, method=label, block=block, coverage=coverage,
                                    common_mask_phase=phase, fp_family_db=scores[0], twofp_family_db=scores[1],
                                    time_map_min_derivative=mind, reconstruction_rejects=reject, frozen_svd_groups=svd_groups))
                for x, y in zip(order, amplitude):
                    spectrum_rows.append(dict(condition=condition, method=label, block=block, carrier_order=x, amplitude=y))
            # These arrays are source data, not a substitute for independent validation.
        np.savez_compressed(folder/'evaluated_events.npz', **saved)
        # Identical event-search freedoms for off-order controls. No score-driven tuning.
        initial_d = contrast_from_energy(energy, motion, cfg)
        for factor in cfg['evaluation']['wrong_event_order_factors']:
            hypothesis = f'wrong_order_{factor:g}'
            print(f'{condition}: {hypothesis} control', flush=True)
            ev = schedule(initial_d, motion, cfg, order_factor=factor)
            dd = contrast_from_energy(energy, motion, cfg, ev)
            ca = candidates(ev, dd, motion, cfg)
            mo = fit_delay_model(ev, ca, cfg['split_seconds']['calibration'], model['smooth'])
            pr = predict_delays(ev, ca, mo, allspans)
            ac = np.isfinite(ca['offsets'][:, :, 0]).sum(axis=1) >= 2
            mask = ac & pr['accepted']
            rw, _ = extract_events(wave, fs, ev, np.zeros((len(ev['t']), 3)))
            te = make_templates(rw, ev, cfg['split_seconds']['calibration'])
            lc = local_xcorr(wave, fs, ev, te, ac)
            for label, shift in [('none', np.zeros_like(lc)), ('local_xcorr', lc), ('constrained', pr['shift'])]:
                ww = rw if label == 'none' else extract_events(wave, fs, ev, shift)[0]
                add_scores(ww, ev, mask, cfg, label, condition, hypothesis, phase_rows)
            save_json(folder/f'{hypothesis}.json', dict(model=mo, origin_rad=ev['origin_rad'],
                                                       accepted_fraction=float(mask.mean()),
                                                       n_events=len(mask)))
        # Encoder scoring starts after vibration-only models and controls are fixed.
        enc = load_npz(folder/'encoder_evaluation_only.npz')
        e = enc['a']
        if len(e) > 100:
            enc_angle = np.interp(motion['time'], e, 2*np.pi*np.arange(len(e))/cfg['encoder_ppr'])
            diff = motion['phase']/cfg['gear_teeth']['ring']-enc_angle
            offset = np.median(diff[interval(motion['time'], cfg['split_seconds']['calibration'])])
            erows = []
            for block, span in enumerate(cfg['split_seconds']['test']):
                use = interval(motion['time'], span)
                erows.append(dict(block=block, angle_rmse_deg=float(np.sqrt(np.mean((diff[use]-offset)**2))*180/np.pi),
                                  edge_rate_per_s=float((len(e)-1)/(e[-1]-e[0]))))
            save_csv(folder/'encoder_evaluation.csv', erows)
        save_csv(out/'metrics.csv', metrics); save_csv(out/'phase_strata.csv', phase_rows)
        save_csv(out/'spectra.csv', spectrum_rows); save_csv(out/'arrival_prediction.csv', arrival_rows)
        save_csv(out/'coverage.csv', coverage_rows)
    # Match stratum and block pairs across target and both controls before Delta S.
    delta_rows = []
    for condition in ['BL', 'PF50']:
        records = [r for r in phase_rows if r['condition'] == condition]
        mapping = {(r['hypothesis'], r['method'], r['branch'], r['angle_bin'], r['channel'], r['block1'], r['block2']):r['phase'] for r in records}
        for method in ['local_xcorr', 'constrained']:
            values = []
            keys = set((r['branch'], r['angle_bin'], r['channel'], r['block1'], r['block2']) for r in records)
            for key in sorted(keys):
                hypotheses = ['target']+[f'wrong_order_{v:g}' for v in cfg['evaluation']['wrong_event_order_factors']]
                if not all((h, m, *key) in mapping for h in hypotheses for m in ['none', method]):
                    continue
                before = mapping[('target', 'none', *key)]-np.mean([mapping[(h, 'none', *key)] for h in hypotheses[1:]])
                after = mapping[('target', method, *key)]-np.mean([mapping[(h, method, *key)] for h in hypotheses[1:]])
                values.append(after-before)
            delta_rows.append(dict(condition=condition, method=method, delta_S=float(np.mean(values)) if values else np.nan,
                                   matched_stratum_pairs=len(values), status='defined' if values else 'insufficient_matched_controls'))
    save_csv(out/'delta_S.csv', delta_rows)
    save_json(out/'evaluation_status.json', dict(complete=True, real_records=2, independent_repeats_per_condition=1,
              status='exploratory_internal_holdout', ringing_phase_enabled=False,
              full_tidal_svd=False, denoising='frozen_local_event_SVD',
              spectra='mean_individual_envelope_on_original_time_axis_before_event_SVD',
              limitation='Phase statistics use common masks but adaptive inference remains; heldout-channel prediction is reported separately.'))
    print('Measured evaluation finished.', flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', action='store_true')
    p.add_argument('--stage', choices=['prepare', 'complete'], default='complete')
    p.add_argument('--resume', type=Path)
    p.add_argument('--reuse-preprocessing', type=Path)
    args = p.parse_args()
    if not args.run:
        p.error('Use --run for explicitly authorized measured processing.')
    if args.resume:
        out = args.resume.resolve()
        cfg = json.loads((out/'runtime_config.json').read_text(encoding='utf-8'))
        if (out/'evaluation_status.json').exists():
            raise FileExistsError('Evaluation already completed; preserve this run.')
    else:
        cfg = json.loads((ROOT/'configs'/'event_alignment_v1_plan.json').read_text(encoding='utf-8'))
        out = ROOT/'results'/('event_alignment_v1_'+datetime.now().strftime('%Y%m%d_%H%M%S'))
        out.mkdir(parents=True, exist_ok=False)
        cfg['run_experiments'] = True
        cfg['status'] = 'implemented_pilot'
        cfg['preprocessing_reused_from'] = str(args.reuse_preprocessing) if args.reuse_preprocessing else None
        cfg['solver_max_iterations'] = 1000
        cfg['implemented_development_selection'] = 'common_band_by_mean_validation_arrival_proxy_plus_5ms_coverage_penalty'
        cfg['implemented_path_penalty_ms_units'] = [1e-5, 1e-5, .03, .03, .48, .48, .3, .3, .3, .3]
        cfg['code_sha256'] = {str(v.relative_to(ROOT)): hashlib.sha256(v.read_bytes()).hexdigest()
                               for v in [ROOT/'src_py'/'event_alignment.py', Path(__file__)]}
        save_json(out/'runtime_config.json', cfg)
    print(f'OUTPUT_DIR={out}', flush=True)
    if not args.resume:
        prepare(out, cfg, reuse=args.reuse_preprocessing)
    if args.stage == 'complete':
        evaluate(out, cfg)


if __name__ == '__main__':
    main()
