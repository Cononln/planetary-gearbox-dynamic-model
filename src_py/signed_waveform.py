"""Signed channel fusion and Fourier diagnostics for a waveform pilot."""
import numpy as np


def signed_summary(tsa, duration, fault_bins=84, carrier_bins=31,
                   mesh_bins=2604, band=(4000., 10000.)):
    """No rectification of the time waveform. FFT magnitudes are spectral amplitudes.

    The comb ratio is exploratory periodic structure, not fault-contact truth or
    physical SNR. Mesh lines and their +/-3 carrier-order neighbourhoods are
    excluded from the fault-spaced comb and its matched background.
    """
    x = np.asarray(tsa, float)
    assert x.ndim == 2 and len(x) == 3 and duration > 0
    centered = x-x.mean(axis=1, keepdims=True)
    fused = centered.mean(axis=0)
    coefficients = np.fft.rfft(fused)/len(fused)
    amplitude = 2*np.abs(coefficients)
    amplitude[0] *= .5
    if len(fused) % 2 == 0:
        amplitude[-1] *= .5
    bins = np.arange(len(amplitude))
    frequency = bins/duration
    in_band = (frequency >= band[0]) & (frequency <= band[1])
    distance_mesh = np.abs(bins-np.rint(bins/mesh_bins)*mesh_bins)
    mesh_mask = in_band & (distance_mesh <= 1)
    centres = bins[in_band & (bins % fault_bins == 0)
                   & (distance_mesh > 3*carrier_bins+10)]
    assert len(centres) > 20
    target_indices = centres[:, None]+np.arange(-1, 2)
    background_indices = centres[:, None]+np.r_[np.arange(-10, -4), np.arange(5, 11)]
    assert background_indices.min() > 0 and background_indices.max() < len(amplitude)
    target_rms = np.sqrt(np.mean(amplitude[target_indices]**2))
    background_rms = np.sqrt(np.mean(amplitude[background_indices]**2))
    channel_energy = np.mean(centered**2)
    fused_energy = np.mean(fused**2)
    metrics = dict(
        signed_fused_rms=float(np.sqrt(fused_energy)),
        mean_channel_rms=float(np.sqrt(channel_energy)),
        channel_coherent_energy_fraction=float(fused_energy/max(channel_energy, 1e-30)),
        fault_spaced_comb_background_db=float(20*np.log10(max(target_rms, 1e-30)/max(background_rms, 1e-30))),
        comb_rms=float(target_rms), matched_background_rms=float(background_rms),
        comb_centres=int(len(centres)),
        mesh_bin_energy_fraction=float(np.sum(amplitude[mesh_mask]**2)/max(np.sum(amplitude[in_band]**2), 1e-30)),
        fp_hz=float(fault_bins/duration), fc_hz=float(carrier_bins/duration),
        fm_hz=float(mesh_bins/duration), fft_grid_hz=float(1/duration),
    )
    return fused, frequency, amplitude, metrics


def peak_preserving_indices(y, first, last, max_bins=1800):
    """Display-only thinning: retain actual min/max samples of every block."""
    width = max(1, int(np.ceil((last-first)/max_bins)))
    selected = []
    for a in range(first, last, width):
        b = min(a+width, last)
        selected.extend([a+int(np.argmin(y[a:b])), a+int(np.argmax(y[a:b]))])
    return np.unique(selected)
