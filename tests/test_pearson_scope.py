from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'runtime'/'phase_python'))
sys.path.insert(0,str(ROOT/'src_py'))
import numpy as np
from pearson_scope import pearson_lag,matrix_tsa


def main():
    rng=np.random.default_rng(7209);x=rng.normal(size=4000)
    moving=3.2*np.roll(x,7)+12
    result,lags,corr=pearson_lag(x,moving,25)
    assert abs(result['lag_samples']-7)<.02 and result['peak']>.999999
    fixed=x[25:-25]
    direct=np.array([np.corrcoef(fixed,moving[25+lag:len(x)-25+lag])[0,1] for lag in lags])
    assert np.max(np.abs(corr-direct))<1e-10
    modified,_,_=pearson_lag(4*x+20,9*moving-30,25)
    assert abs(modified['lag_samples']-result['lag_samples'])<1e-8
    z=np.ones(4000);constant,_,_=pearson_lag(z,z,25)
    assert not constant['valid']
    y=rng.normal(size=(3,800))+rng.normal(size=(3,1))
    centered=y-y.mean(axis=1,keepdims=True)
    u,s,v=np.linalg.svd(centered,full_matrices=False)
    for rank in [1,2,3]:
        got,saved=matrix_tsa(y,rank)
        expected=((u[:,:rank]*s[:rank])@v[:rank]+y.mean(axis=1,keepdims=True)).mean(axis=0)
        assert np.allclose(got,expected)
    assert np.allclose(matrix_tsa(y,3)[0],matrix_tsa(y,0)[0])
    # One lag per segment follows known changing channel delay.
    delays=[-9,0,11]
    parts=[rng.normal(size=6000) for _ in delays]
    other=[np.roll(a,lag) for a,lag in zip(parts,delays)]
    estimates=[pearson_lag(a,b,20)[0]['lag_samples'] for a,b in zip(parts,other)]
    assert np.max(np.abs(np.array(estimates)-delays))<.03
    whole=pearson_lag(np.concatenate(parts),np.concatenate(other),20)[0]['lag_samples']
    print('PASS: exact Pearson formula, lag sign, gain/DC invariance, zero-variance rejection, full SVD equivalence, rank-3 identity')
    print('Known segment lags',delays,'estimated',estimates,'one whole-record lag',whole)


if __name__=='__main__':main()
