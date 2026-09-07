"""Fixed-support Pearson lag search and full-tidal matrix averaging."""
import numpy as np
from scipy import signal


def pearson_lag(reference, moving, max_lag):
    """Positive lag means corrected moving[t] = moving[t + lag].

    Reference support is fixed for every candidate lag. Correlation is the
    signed Pearson coefficient, not an absolute value or a raw dot product.
    """
    a=np.asarray(reference,float);b=np.asarray(moving,float)
    assert len(a)==len(b) and len(a)>2*max_lag+8 and max_lag>=1
    fixed=a[max_lag:-max_lag];fixed=fixed-fixed.mean()
    n=len(fixed)
    sums=np.r_[0.,np.cumsum(b)];squares=np.r_[0.,np.cumsum(b*b)]
    variance=np.maximum(squares[n:]-squares[:-n]-(sums[n:]-sums[:-n])**2/n,0)
    denominator=np.sqrt(np.sum(fixed*fixed)*variance)
    numerator=signal.correlate(b,fixed,mode='valid',method='fft')
    corr=np.full(len(numerator),np.nan)
    np.divide(numerator,denominator,out=corr,where=denominator>1e-18)
    lags=np.arange(-max_lag,max_lag+1)
    if not np.isfinite(corr).any():
        return dict(lag_samples=0.,peak=np.nan,zero=np.nan,boundary=False,valid=False),lags,corr
    ix=int(np.nanargmax(corr));sub=0.
    if 0<ix<len(corr)-1 and np.isfinite(corr[ix-1:ix+2]).all():
        den=corr[ix-1]-2*corr[ix]+corr[ix+1]
        if den<0:sub=float(np.clip(.5*(corr[ix-1]-corr[ix+1])/den,-.5,.5))
    return dict(lag_samples=float(lags[ix]+sub),peak=float(corr[ix]),zero=float(corr[max_lag]),
                boundary=ix in [0,len(corr)-1],valid=True),lags,corr


def pearson(a,b):
    x=np.asarray(a,float);y=np.asarray(b,float)
    x=x-x.mean();y=y-y.mean()
    den=np.sqrt(np.sum(x*x)*np.sum(y*y))
    return float(np.sum(x*y)/den) if den>1e-20 else np.nan


def matrix_tsa(rows,rank):
    """Exact low-rank row-space SVD projection followed by row average.

    Gram eigendecomposition avoids a very wide full SVD. This is an
    in-sample denoiser; the output must not be labelled held-out validation.
    """
    x=np.asarray(rows,float);dc=x.mean(axis=1,keepdims=True);z=x-dc
    if rank==0:
        return x.mean(axis=0),np.array([])
    assert 1<=rank<=len(x)
    eigen,u=np.linalg.eigh(z@z.T)
    order=np.argsort(eigen)[::-1];eigen=np.maximum(eigen[order],0);u=u[:,order[:rank]]
    weights=np.ones(len(x))/len(x)@u@u.T
    mean=weights@z+dc.mean()
    return mean,np.sqrt(eigen)
