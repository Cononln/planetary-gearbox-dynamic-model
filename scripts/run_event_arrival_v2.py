"""Measured envelope-arrival pilot with development-frozen candidate constraints."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path
import os
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'runtime'/'phase_python'))
sys.path.insert(0, str(ROOT/'src_py'))
os.environ.setdefault('OMP_NUM_THREADS', '2')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
import numpy as np
from event_alignment import (interval, schedule, contrast_from_energy, candidates, cross_block_metrics,
                              extract_events, reconstruct_block, family_spectrum)
from event_arrival import (envelope_response, onset_markers, arrival_templates, envelope_candidates,
                           local_result, fit_relative_path, joint_result, other_channel_baseline,
                           cross_channel_marker_error)


def load_npz(path):
    with np.load(path) as archive:
        return {key: archive[key] for key in archive.files}


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, default=jsonable)+'\n', encoding='utf-8')


def save_csv(path, rows):
    if not rows:
        return
    with path.open('w', newline='', encoding='utf-8-sig') as stream:
        w = csv.DictWriter(stream, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def heldout(events, cand, marker, model, spans, strength, condition, phase):
    rows = []
    for j in range(3):
        joint = joint_result(events, cand, model, spans, strength, excluded=j)
        baseline, valid = other_channel_baseline(events, cand, model, j)
        for block, span in enumerate(spans):
            eligible = interval(events['t'], span)
            observed = np.isfinite(marker['offset'][:, j])
            common = eligible & observed & valid & joint['accepted']
            for method, prediction, own in [('none', np.zeros(len(eligible)), valid),
                    ('envelope_other_median', baseline, valid),
                    ('joint_other_prediction', joint['predicted'][:, j], joint['accepted'])]:
                error = np.abs(prediction[common]-marker['offset'][common, j])*1000
                own_mask = eligible & observed & own
                rows.append(dict(condition=condition, phase=phase, strength=strength, channel=j+1, block=block,
                                 method=method, mae_ms=float(error.mean()) if len(error) else np.nan,
                                 count=int(common.sum()), own_coverage=float(own_mask.sum()/max(1, eligible.sum())),
                                 common_coverage=float(common.sum()/max(1, eligible.sum()))))
    return rows


def prepare_condition(source, condition, cfg, out):
    folder = source/condition
    fs = cfg['working_rates_hz']['impulse']
    events = load_npz(folder/'band_2'/'events.npz')
    support = np.isfinite(load_npz(folder/'band_2'/'candidates.npz')['offsets'][:, :, 0])
    wave = np.load(folder/'band_2'/'waveform.npy', mmap_mode='r')
    envelope = envelope_response(wave, fs)
    marker = onset_markers(envelope, fs, events)
    templates, audit = arrival_templates(envelope, fs, events, marker, cfg['split_seconds']['calibration'])
    cand = envelope_candidates(envelope, fs, events, templates, support)
    model = fit_relative_path(events, cand, cfg['split_seconds']['calibration'])
    target = out/condition; target.mkdir()
    np.save(target/'envelope.npy', envelope)
    np.savez_compressed(target/'events.npz', **events)
    np.savez_compressed(target/'markers.npz', **marker)
    np.savez_compressed(target/'envelope_candidates.npz', **cand)
    np.savez_compressed(target/'templates.npz', **{f'{c}_{b}_{j}': value for (c,b,j), value in templates.items()})
    save_json(target/'template_audit.json', audit); save_json(target/'path_model.json', model)
    results = []
    for strength in [.2, .5, 1.]:
        results += heldout(events, cand, marker, model, [cfg['split_seconds']['validation']], strength, condition, 'validation')
    print(f'{condition}: prepared, {len(events["t"])} schedules; {model["calibration_supported_events"]} calibration events', flush=True)
    return results


def evaluate_condition(source, condition, cfg, out, strength):
    folder = out/condition; fs = cfg['working_rates_hz']['impulse']
    events = load_npz(folder/'events.npz'); marker = load_npz(folder/'markers.npz')
    cand = load_npz(folder/'envelope_candidates.npz')
    model = json.loads((folder/'path_model.json').read_text(encoding='utf-8'))
    env = np.load(folder/'envelope.npy', mmap_mode='r')
    wave = np.load(source/condition/'band_2'/'waveform.npy', mmap_mode='r')
    motion = load_npz(source/condition/'motion.npz')
    test = cfg['split_seconds']['test']
    spans = [cfg['split_seconds']['calibration'], cfg['split_seconds']['validation']]+test
    local = local_result(cand)
    joint = joint_result(events, cand, model, spans, strength)
    with np.load(source/condition/'evaluated_events.npz') as archive:
        v1 = dict(shift=archive['constrained_shift'], accepted=archive['constrained_mask'])
        rawcorr = dict(shift=archive['local_xcorr_shift'], accepted=archive['local_xcorr_mask'])
    methods = dict(none=dict(shift=np.zeros_like(local['shift']), accepted=local['accepted']),
                   waveform_xcorr=rawcorr, v1_constraint=v1, envelope_xcorr=local, joint_envelope=joint)
    common = local['accepted'] & joint['accepted'] & rawcorr['accepted'] & v1['accepted']
    main_common = local['accepted'] & joint['accepted']
    save_json(folder/'coverage_masks.json', dict(main='local envelope & joint envelope', all_five='intersection including V1; older V1 rejection may reduce support'))
    low = np.load(source/condition/'band_0'/'waveform.npy', mmap_mode='r')
    lowenv = envelope_response(low, fs)
    lowmarker = onset_markers(lowenv, fs, events)
    np.savez_compressed(folder/'cross_band_markers.npz', **lowmarker)
    del lowenv
    arrival = heldout(events, cand, marker, model, test, strength, condition, 'test')
    metrics, strata, spectra = [], [], []
    saved = dict(t=events['t'], branch=events['branch'], psi=events['psi'], common_mask=common,
                 main_mask=main_common, marker_offsets=marker['offset'], cross_band_marker_offsets=lowmarker['offset'],
                 joint_assignment=joint['assignment'], path=joint['path'], common=joint['common'])
    for name, result in methods.items():
        print(f'{condition}: evaluating {name}', flush=True)
        windows, u = extract_events(env, fs, events, result['shift'], core=(-.012, .024))
        rows = cross_block_metrics(windows, events, common, test)
        primary_rows = cross_block_metrics(windows, events, main_common, test) if name in ['none','envelope_xcorr','joint_envelope'] else []
        strata += [dict(condition=condition, hypothesis='target', method=name, mask='all_five', **r) for r in rows]
        strata += [dict(condition=condition, hypothesis='target', method=name, mask='main_three', **r) for r in primary_rows]
        score = float(np.mean([r['phase'] for r in rows])) if rows else np.nan
        main_score = float(np.mean([r['phase'] for r in primary_rows])) if primary_rows else np.nan
        saved[name] = windows; saved[name+'_shift'] = result['shift']; saved[name+'_mask'] = result['accepted']; saved['u'] = u
        for block, span in enumerate(test):
            testmask = interval(events['t'], span)
            value, pairs = cross_channel_marker_error(lowmarker, result['shift'], common & testmask)
            inband, in_pairs = cross_channel_marker_error(marker, result['shift'], common & testmask)
            signal_block, derivative, rejected = reconstruct_block(wave, fs, events, result['shift'], result['accepted'], span)
            order, amplitude, family = family_spectrum(signal_block, fs, span, motion)
            metrics.append(dict(condition=condition, method=name, block=block,
                           envelope_consistency=score, main_envelope_consistency=main_score,
                           cross_band_arrival_spread_ms=value, cross_band_pairs=pairs,
                           inband_arrival_spread_ms=inband, inband_pairs=in_pairs,
                           coverage=float(result['accepted'][testmask].mean()),
                           common_events=int((common & testmask).sum()), main_events=int((main_common & testmask).sum()),
                           family_db=family[0], twice_family_db=family[1], time_map_min_derivative=derivative,
                           reconstruction_rejected=rejected))
            spectra += [dict(condition=condition, method=name, block=block, carrier_order=x, amplitude=y) for x,y in zip(order, amplitude)]
    np.savez_compressed(folder/'evaluated_envelopes.npz', **saved)
    energy = np.load(source/condition/'band_2'/'energy.npy', mmap_mode='r')
    d0 = contrast_from_energy(energy, motion, cfg)
    for factor in [.85, 1.15]:
        print(f'{condition}: wrong-order {factor}', flush=True)
        ev = schedule(d0, motion, cfg, order_factor=factor)
        dd = contrast_from_energy(energy, motion, cfg, ev)
        evidence_cand = candidates(ev, dd, motion, cfg)
        support = np.isfinite(evidence_cand['offsets'][:, :, 0])
        mark = onset_markers(env, fs, ev)
        templates, _ = arrival_templates(env, fs, ev, mark, cfg['split_seconds']['calibration'])
        ca = envelope_candidates(env, fs, ev, templates, support)
        mo = fit_relative_path(ev, ca, cfg['split_seconds']['calibration'])
        lc = local_result(ca); pr = joint_result(ev, ca, mo, spans, strength)
        mask = lc['accepted'] & pr['accepted']
        for name, result in [('none',dict(shift=np.zeros_like(lc['shift']))),('envelope_xcorr',lc),('joint_envelope',pr)]:
            ww, _ = extract_events(env, fs, ev, result['shift'], core=(-.012,.024))
            rr = cross_block_metrics(ww, ev, mask, test)
            strata += [dict(condition=condition, hypothesis=f'wrong_{factor}', method=name, mask='main_three', **r) for r in rr]
        np.savez_compressed(folder/f'wrong_{factor}.npz', t=ev['t'], branch=ev['branch'], psi=ev['psi'],
                            envelope_shift=lc['shift'], joint_shift=pr['shift'], accepted=mask)
        save_json(folder/f'wrong_{factor}_model.json', mo)
    return metrics, strata, spectra, arrival


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--run', action='store_true')
    parser.add_argument('--source', type=Path, default=ROOT/'results'/'event_alignment_v1_20260906_053745')
    args = parser.parse_args()
    if not args.run:
        parser.error('Use --run to execute the measured pilot.')
    source = args.source.resolve()
    cfg = json.loads((source/'runtime_config.json').read_text(encoding='utf-8'))
    out = ROOT/'results'/('event_arrival_v2_'+datetime.now().strftime('%Y%m%d_%H%M%S'))
    out.mkdir()
    status = dict(version='event_arrival_v2', complete=False, source=str(source),
                  inherited_motion_and_band=True, primary='heldout event-onset proxy plus known-arrival controls',
                  amplitude_inverse_compensation=False, ringing_phase_alignment=False,
                  envelope_smoothing_sigma_s=.00012, onset_fraction=.2, onset_search_s=.003,
                  ncc_minimum=.25, candidates_per_channel=4, lag_search_s=.008, candidate_spacing_s=.0008,
                  path_strength_candidates=[.2,.5,1.], dp_common_grid_ms=[-8,8,.2],
                  dp_transition_strength=.08, timing_penalty_scale_ms=.7, missing_cost=2.5,
                  event_display_and_score_core_s=[-.012,.024], template_core_s=[-.004,.016],
                  split_seconds=cfg['split_seconds'], independent_records_per_condition=1,
                  evidence_status='exploratory_internal_holdout_on_historically_inspected_records',
                  code_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in [ROOT/'src_py'/'event_arrival.py', ROOT/'src_py'/'event_alignment.py', Path(__file__)]})
    save_json(out/'run_manifest.json', status)
    print(f'OUTPUT_DIR={out}', flush=True)
    development = []
    for condition in ['BL','PF50']:
        development += prepare_condition(source, condition, cfg, out)
    save_csv(out/'validation_prediction.csv', development)
    losses = []
    for strength in [.2,.5,1.]:
        rr = [r for r in development if r['strength']==strength and r['method']=='joint_other_prediction']
        vals = [r['mae_ms']+5*(1-r['own_coverage']) if np.isfinite(r['mae_ms']) else 50. for r in rr]
        losses.append(dict(strength=strength, validation_loss_ms=float(np.mean(vals))))
    selected = min(losses, key=lambda r:r['validation_loss_ms'])['strength']
    save_json(out/'selection.json', dict(strength=selected, candidates=losses,
                                       criterion='mean held-out onset proxy MAE + 5 ms*(1-coverage), both states, validation only'))
    print(f'FROZEN_STRENGTH={selected}', flush=True)
    allrows = [[],[],[],[]]
    for condition in ['BL','PF50']:
        result = evaluate_condition(source, condition, cfg, out, selected)
        for dest, rows in zip(allrows,result):
            dest.extend(rows)
        for name, rows in zip(['metrics','envelope_strata','spectra','heldout_prediction'],allrows):
            save_csv(out/(name+'.csv'),rows)
    ds = []
    for condition in ['BL','PF50']:
        rr = [r for r in allrows[1] if r['condition']==condition and r['mask']=='main_three']
        keys = ['branch','angle_bin','channel','block1','block2']
        mapping = {(r['hypothesis'],r['method'],*(r[k] for k in keys)):r['phase'] for r in rr}
        strata = sorted(set(tuple(r[k] for k in keys) for r in rr))
        for method in ['envelope_xcorr','joint_envelope']:
            values=[]
            for key in strata:
                if not all((h,m,*key) in mapping for h in ['target','wrong_0.85','wrong_1.15'] for m in ['none',method]):
                    continue
                target = mapping[('target',method,*key)]-mapping[('target','none',*key)]
                control = np.mean([mapping[(h,method,*key)]-mapping[(h,'none',*key)] for h in ['wrong_0.85','wrong_1.15']])
                values.append(target-control)
            ds.append(dict(condition=condition,method=method,delta_envelope_gain=float(np.mean(values)) if values else np.nan,
                           matched_stratum_pairs=len(values)))
    save_csv(out/'wrong_order_gain.csv',ds)
    status['complete']=True; status['selected_strength']=selected
    save_json(out/'run_manifest.json',status)
    print('Completed measured V2 pilot.',flush=True)


if __name__=='__main__':
    main()
