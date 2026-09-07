from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'runtime'/'phase_python'))
sys.path.insert(0, str(ROOT/'src_py'))
import numpy as np
from scipy import signal
from pearson_scope import pearson_lag, matrix_tsa
from signed_waveform import signed_summary, peak_preserving_indices


def forbidden(*args, **kwargs):
    raise AssertionError('Unexpected analytic-envelope operation')


def main():
    signal.hilbert = forbidden
    rng = np.random.default_rng(51200)
    n = 65536
    duration = 2.
    t = np.arange(n)/n*duration
    sine = 2*np.sin(2*np.pi*6000*t)
    _, f, amp, metric = signed_summary(np.stack([sine]*3), duration)
    assert np.isclose(amp[np.argmin(abs(f-6000))], 2)
    assert np.isclose(metric['channel_coherent_energy_fraction'], 1)
    cancelled, _, _, m = signed_summary(np.stack([sine, -sine, np.zeros(n)]), duration)
    assert np.max(abs(cancelled)) < 1e-12 and m['signed_fused_rms'] < 1e-12
    x = rng.normal(size=6000)
    result, _, _ = pearson_lag(x, np.roll(x, 9), 20)
    assert abs(result['lag_samples']-9) < .02
    rows = rng.normal(size=(3, 6000))
    assert np.allclose(matrix_tsa(rows, 3)[0], rows.mean(axis=0))
    extreme = np.zeros(200)
    extreme[42], extreme[43] = 30, -29
    keep = peak_preserving_indices(extreme, 0, 200, 10)
    assert 42 in keep and 43 in keep
    # Synchronous interference survives averaging; independent noise attenuates.
    target = np.exp(-((t-.3)/.003)**2)*np.cos(2*np.pi*6000*(t-.3))
    synchronous_mesh = .8*np.sin(2*np.pi*7000*t)
    noisy = target+synchronous_mesh+rng.normal(scale=.4, size=(200, n))
    residual = noisy.mean(axis=0)-target-synchronous_mesh
    assert np.std(residual) < .04
    print('PASS: signed cancellation, FFT amplitude, gain-normalized lag, SVD identity, '
          'peak-preserving display, synchronous-interference survival, analytic-envelope guard')


if __name__ == '__main__':
    main()
