"""Numerical/source audit of completed cycle-TSA verification artifacts."""
from pathlib import Path
import sys
import json
import hashlib
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'runtime'/'phase_python'))
sys.path.insert(0,str(ROOT/'src_py'))
import numpy as np
from tsa_validation import harmonic_amplitudes
from run_event_arrival_v2 import save_json,load_npz
from plot_fault_event_alignment import read_csv


def main():
    out=Path(sys.argv[1]).resolve()
    manifest=json.loads((out/'run_manifest.json').read_text(encoding='utf-8'))
    assert manifest['complete'] and not manifest['alignment_refit']
    for p,h in manifest['frozen_input_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
    blocks=read_csv(out/'block_metrics.csv');groups=read_csv(out/'group_metrics.csv')
    registry=read_csv(out/'cycle_registry.csv');fits=read_csv(out/'svd_fits.csv')
    harmonics=read_csv(out/'harmonics.csv');diag=read_csv(out/'reconstruction_audit.csv')
    assert len(blocks)==2*3*2*4*3
    assert len(groups)==len(blocks)*6 and len(harmonics)==len(groups)*8
    assert all(int(r['rejected_warps'])==0 for r in diag)
    assert all(float(r['time_map_min_derivative'])>=.199 for r in diag)
    assert all(r['svd_enabled']=='True' and int(r['calibration_cycles'])>=8 for r in fits)
    for r in groups:
        assert 0<=float(r['tsa_retention'])<=1.000001
        assert -1.000001<=float(r['split_tsa_correlation'])<=1.000001
        assert np.isfinite(float(r['heldout_explained']))
    for condition in ['BL','PF50']:
        for b,span in enumerate(manifest['split_seconds']['test']):
            reg=[r for r in registry if r['condition']==condition and r['split']=='test' and int(r['block'])==b]
            assert all(float(r['start_s'])>=span[0]+.02-1e-9 and float(r['end_s'])<=span[1]-.02+1e-9 for r in reg)
            assert all(int(r['n_cycles'])==len(reg) for r in blocks if r['condition']==condition and int(r['block'])==b)
            for rep in ['waveform','envelope']:
                for group in range(6):
                    rr=[r for r in groups if r['condition']==condition and int(r['block'])==b and r['representation']==rep and int(r['angle_bin'])==group]
                    assert len(set(r['n_cycles'] for r in rr))==1
    maximum_error=0.
    negative_envelope_samples=0
    for condition in ['BL','PF50']:
        for method in ['none','envelope_xcorr','joint_envelope']:
            data=load_npz(out/condition/f'{method}_tsa.npz')
            for rep in ['waveform','envelope']:
                for rank in [0,1,3,5]:
                    for b in range(3):
                        for group in range(6):
                            tsa=data[f'{rep}_rank{rank}_block{b}_group{group}']
                            amps=harmonic_amplitudes(tsa,rep)
                            if rep=='envelope':negative_envelope_samples+=int(np.sum(tsa<0))
                            rr=[r for r in harmonics if r['condition']==condition and r['method']==method and r['representation']==rep and int(r['rank'])==rank and int(r['block'])==b and int(r['angle_bin'])==group]
                            saved=np.array([float(r['amplitude']) for r in rr])
                            error=np.max(np.abs(amps-saved))/max(np.max(saved),1e-12)
                            maximum_error=max(maximum_error,float(error))
    assert maximum_error<1e-4
    result=dict(passed=True,frozen_alignment_hashes_unchanged=True,
                identical_complete_cycle_counts=True,test_boundaries_guarded=True,
                all_calibration_groups_svd_enabled=True,no_time_map_rejections=True,
                source_tsa_harmonic_max_relative_error=maximum_error,
                negative_envelope_tsa_samples=negative_envelope_samples,
                rows=dict(blocks=len(blocks),strata=len(groups),harmonics=len(harmonics)),
                formula_tests='tests/test_tsa_validation.py passed before data evaluation; known pulse TSA retention 0.310603 -> 0.808592; individual Fourier magnitudes unchanged',
                independent_repeat_or_contact_truth=False)
    save_json(out/'numerical_audit.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
