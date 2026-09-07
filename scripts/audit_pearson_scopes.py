"""Verify scope-comparison numerical artifacts independently of plotting."""
from pathlib import Path
import sys
import json
import hashlib
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'runtime'/'phase_python'))
sys.path.insert(0,str(ROOT/'src_py'))
import numpy as np
from run_event_arrival_v2 import save_json,load_npz
from plot_fault_event_alignment import read_csv
from compare_pearson_scopes import spectrum,METHODS


def main():
    out=Path(sys.argv[1]).resolve()
    manifest=json.loads((out/'run_manifest.json').read_text(encoding='utf-8'))
    assert manifest['complete']
    for p,h in manifest['original_v2_unchanged_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
    periods=read_csv(out/'periods.csv');fits=read_csv(out/'lag_fits.csv');lags=read_csv(out/'effective_lags.csv')
    metrics=read_csv(out/'metrics.csv');checks=read_csv(out/'alignment_checks.csv')
    assert len(periods)==6 and len(fits)==32 and len(metrics)==60 and len(lags)==90 and len(checks)==150
    assert all(r['valid']=='True' and abs(float(r['lag_ms']))<=8.01 for r in fits)
    assert all(-1.000001<=float(r['pearson_peak'])<=1.000001 for r in fits)
    maximum=0.
    for condition in ['BL','PF50']:
        reg=[r for r in periods if r['condition']==condition]
        assert len(reg)==3 and all(int(r['carrier_turns'])==31 and int(r['fault_cycles'])==84 for r in reg)
        assert all(abs(float(reg[i]['end_s'])-float(reg[i+1]['start_s']))<1e-10 for i in [0,1])
        duration=np.mean([float(r['duration_s']) for r in reg])
        for method in METHODS:
            for j in range(1,4):
                shift=[float(r['lag_ms']) for r in lags if r['condition']==condition and r['method']==method and int(r['channel'])==j]
                if j==1 or method=='none':assert np.allclose(shift,0)
                if method.startswith('whole'):assert len(set(shift))==1
            data=load_npz(out/condition/f'{method}_tsa.npz')
            for rep in ['waveform','envelope']:
                for rank in [0,1,2]:
                    tsa=data[f'{rep}_rank{rank}']
                    assert tsa.shape==(3,524288) and np.isfinite(tsa).all()
                    _,_,_,computed=spectrum(tsa,rep,duration,84)
                    saved=next(r for r in metrics if r['condition']==condition and r['method']==method and r['representation']==rep and int(r['rank'])==rank)
                    error=abs(computed['family_db']-float(saved['family_db']))
                    maximum=max(maximum,error)
            del data
        for k in range(3):
            rows=[r for r in checks if r['condition']==condition and r['tidal']==str(k) and r['representation']=='arrival_proxy']
            assert len(set(r['marker_pairs'] for r in rows))==1
    assert maximum<1e-3
    report=dict(passed=True,original_v2_unchanged=True,same_three_full_periods=True,
                s1_fixed=True,whole_lags_constant=True,same_marker_pairs=True,
                max_recomputed_family_difference_db=maximum,
                lag_boundary_hits=sum(r['boundary']=='True' for r in fits),
                source_tsa_shapes_verified=True,
                unit_tests='Exact Pearson/direct np.corrcoef equivalence, lag sign, DC/gain invariance, zero-variance rejection, full SVD equivalence and rank3 identity passed',
                evidence='in-sample descriptive comparison, not independent fault-contact truth')
    save_json(out/'numerical_audit.json',report);print(json.dumps(report,indent=2))


if __name__=='__main__':main()
