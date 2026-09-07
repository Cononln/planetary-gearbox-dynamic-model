from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'runtime'/'phase_python'))
sys.path.insert(0, str(ROOT/'src_py'))
import numpy as np
from tidal_waveform import recurrence, angular_rows, split_row_correlation
from signed_waveform import signed_summary


def main():
    sun = recurrence('sun')
    planet = recurrence('planet')
    assert sun == dict(carrier_turns=1, mesh_cycles=84, fault_cycles=12)
    assert planet == dict(carrier_turns=31, mesh_cycles=2604, fault_cycles=84)
    t = np.arange(0, 60, .001)
    edges, q = angular_rows(dict(time=t, phase=2*np.pi*168*t), 1, 16384)
    assert q.shape == (112, 16384)
    assert np.allclose(np.diff(edges), .5)
    assert np.allclose(q[:, 0], edges[:-1])
    assert np.all(q[:, -1] < edges[1:])
    z = np.sin(2*np.pi*6000*np.arange(16384)/32768)
    assert np.isclose(split_row_correlation(np.tile(z, (20, 1))), 1)
    _, f, a, m = signed_summary(np.tile(z, (3, 1)), .5,
                                fault_bins=12, carrier_bins=1, mesh_bins=84)
    assert np.isclose(a[np.argmin(abs(f-6000))], 1)
    assert m['fp_hz'] == 24 and m['fm_hz'] == 168 and m['fft_grid_hz'] == 2
    print('PASS: sun/planet recurrence, sun row count and duration, endpoint convention, signed split repeatability, sun spectral bins')


if __name__ == '__main__':
    main()
