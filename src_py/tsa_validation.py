"""Frozen-basis cycle TSA and cross-validated diagnostics, not new estimators."""
from __future__ import annotations
import numpy as np
from scipy import signal


def fit_basis(rows, max_rank=5):
    x = np.asarray(rows, float)
    x = x - x.mean(axis=-1, keepdims=True)
    _, singular, right = np.linalg.svd(x, full_matrices=False)
    return right[:max_rank], singular


def project(rows, basis, rank):
    x = np.asarray(rows, float)
    if rank == 0 or basis is None:
        return x.copy()
    dc = x.mean(axis=-1, keepdims=True)
    v = basis[:rank]
    return ((x-dc) @ v.T) @ v + dc


def envelope_of_tsa(tsa, representation):
    if representation == 'waveform':
        return np.abs(signal.hilbert(tsa, axis=-1))
    if representation == 'envelope':
        return np.asarray(tsa)
    raise ValueError(representation)


def harmonic_amplitudes(tsa, representation, n_harmonics=8):
    """Mean channel envelopes, then Fourier series on a complete cycle."""
    env = envelope_of_tsa(tsa, representation).mean(axis=0)
    return 2*np.abs(np.fft.rfft(env-env.mean()))[1:n_harmonics+1]/len(env)


def tsa_metrics(original, processed):
    """Input tensors: events x channels x cycle samples, same row order.

    The test target in prediction remains original even with SVD enabled.
    Scalar row DC is removed throughout the energy diagnostics.
    """
    x = np.asarray(original, float)
    y = np.asarray(processed, float)
    assert x.shape == y.shape and x.ndim == 3 and len(x) >= 4
    x = x-x.mean(axis=-1, keepdims=True)
    y = y-y.mean(axis=-1, keepdims=True)
    input_energy = float(np.mean(x*x))
    output_energy = float(np.mean(y.mean(axis=0)**2))
    n = (len(x)//2)*2
    xa, xb = x[:n:2], x[1:n:2]
    ya, yb = y[:n:2].mean(axis=0), y[1:n:2].mean(axis=0)
    cross_prediction = .5*(np.mean(2*xb.mean(axis=0)*ya-ya*ya)
                             + np.mean(2*xa.mean(axis=0)*yb-yb*yb))
    paired_input = .5*(np.mean(xa*xa)+np.mean(xb*xb))
    norm = np.sqrt(np.sum(ya*ya)*np.sum(yb*yb))
    return dict(n_cycles=len(x), paired_cycles=n,
                input_ac_rms=np.sqrt(input_energy), tsa_ac_rms=np.sqrt(output_energy),
                tsa_retention=output_energy/max(input_energy, 1e-30),
                heldout_explained=float(cross_prediction/max(paired_input, 1e-30)),
                split_tsa_correlation=float(np.sum(ya*yb)/max(norm, 1e-30)))
