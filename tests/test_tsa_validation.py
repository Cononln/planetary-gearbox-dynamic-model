"""Known-signal checks for downstream TSA; no empirical superiority assumptions."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'runtime'/'phase_python'))
sys.path.insert(0,str(ROOT/'src_py'))
sys.path.insert(0,str(ROOT/'scripts'))
import numpy as np
from tsa_validation import fit_basis, project, tsa_metrics, harmonic_amplitudes


def main():
    rng=np.random.default_rng(7106)
    q=np.arange(2048)/2048
    pulse=np.exp(-.5*((q-.4)/.012)**2)
    clean=2+pulse
    shifts=np.linspace(-90,90,60).astype(int)
    rows=np.stack([np.roll(clean,s) for s in shifts])
    noise=rng.normal(0,.07,rows.shape)
    raw=np.repeat((rows+noise)[:,None,:],3,axis=1)
    aligned=np.stack([np.roll(raw[i],-s,axis=-1) for i,s in enumerate(shifts)])
    a=tsa_metrics(raw,raw);b=tsa_metrics(aligned,aligned)
    assert b['tsa_retention'] > a['tsa_retention']+.15
    assert b['heldout_explained'] > a['heldout_explained']+.15
    # Whole-row translations preserve individual Fourier amplitudes.
    assert np.allclose(np.abs(np.fft.rfft(raw)),np.abs(np.fft.rfft(aligned)),atol=1e-11)
    signs=(-1.)**np.arange(60)
    ringing=pulse*np.cos(2*np.pi*130*q)
    w=signs[:,None,None]*np.broadcast_to(ringing,(60,3,len(q)))
    assert np.max(np.abs(w.mean(axis=0)))<1e-12
    assert tsa_metrics(w,w)['tsa_retention']<1e-20
    env=np.broadcast_to(clean,(60,3,len(q))).copy()
    assert tsa_metrics(env,env)['tsa_retention']>.999999
    cal=rng.normal(size=(12,len(q)))+pulse
    v,_=fit_basis(cal)
    immutable=v.copy()
    test=rng.normal(size=(8,len(q)))+pulse
    out=project(test,v,3)
    assert np.allclose(out.mean(axis=0),project(test.mean(axis=0,keepdims=True),v,3)[0])
    assert np.allclose(out.mean(axis=-1),test.mean(axis=-1))
    changed=test.copy();changed[4:]*=100
    assert np.allclose(project(changed,v,3)[:4],out[:4])
    assert np.array_equal(v,immutable)
    assert np.allclose(project(test,v,0),test)
    sine=5+2.3*np.cos(2*np.pi*3*q+.7)
    amps=harmonic_amplitudes(np.stack([sine]*3),'envelope')
    assert abs(amps[2]-2.3)<1e-12
    assert np.max(np.delete(amps,2))<1e-12
    orthogonal=np.eye(8,2048)[:,None,:].repeat(3,axis=1)
    assert tsa_metrics(orthogonal,orthogonal)['heldout_explained']<0
    from validate_tsa_v2 import cycle_registry
    t=np.linspace(0,10,10001);fp=5.0;zp=31
    event_t=np.arange(.2,9.9,.1)
    events=dict(t=event_t,branch=np.arange(len(event_t))%2,psi=2*np.pi*event_t)
    motion=dict(time=t,phase=2*np.pi*fp*zp*t)
    accepted=np.ones(len(event_t),bool)
    registry=cycle_registry(events,accepted,motion,[1,3],zp)
    assert len(registry)>5
    assert all(abs(r['end_s']-r['start_s']-1/fp)<1e-10 for r in registry)
    assert all(r['start_s']>=1.02 and r['end_s']<=2.98 for r in registry)
    bad_index=registry[2]['event_index']+1
    accepted[bad_index]=False
    rejected=cycle_registry(events,accepted,motion,[1,3],zp)
    assert len(rejected)==len(registry)-1
    print('PASS: known shifts, unchanged individual spectra, signed cancellation, envelope TSA, Fourier scaling, frozen SVD, linearity, negative heldout control')
    print('PASS: complete-cycle period, split guards, shared two-branch support')
    print('Known pulse TSA retention:',a['tsa_retention'],'->',b['tsa_retention'])
    print('Known pulse heldout explained:',a['heldout_explained'],'->',b['heldout_explained'])


if __name__=='__main__':main()
