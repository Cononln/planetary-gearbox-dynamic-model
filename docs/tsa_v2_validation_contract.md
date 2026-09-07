# V2 downstream TSA verification

Question: do frozen envelope-correlation or joint event shifts improve synchronous averaging, and does a fixed SVD denoiser change that comparison?

Backend: Python. Figure archetype: quantitative grids. Output: 183 mm figures, editable SVG/PDF, PNG and 600 dpi TIFF, source tables and run manifest. Provisional conclusion is determined by the results, not assumed in advance.

## Fixed protocol

- Reuse existing BL and PF50 recordings, frozen motion estimate, event origins, shifts and the 4–10 kHz response. No new alignment fitting, rank selection, encoder use or amplitude compensation.
- Compare none, envelope_xcorr, joint_envelope on a common acceptance mask. Use the same complete cycles in each arm. Require both event branches within each cycle to be supported.
- One matrix row spans one full same-side planet-tooth recurrence, mesh phase / 31 changing by 2 pi. Use phase interval [-0.25, 0.75) relative to branch A so both A and B contacts lie inside the row. Keep only rows contained in the original calibration/test partitions with an additional 20 ms guard. Resample to 12,288 angular samples per cycle.
- Preserve six 60-degree carrier-angle strata and three separate sensor channels. This is fault-cycle, path-conditioned TSA, not a full 31-carrier-turn tidal matrix. Each 10 s test partition is shorter than the approximately 15.6 s full tidal recurrence.
- Two explicitly separate representations: signed band vibration (SVD/TSA then Hilbert envelope) and per-cycle envelope (Hilbert envelope before SVD/TSA). Neither is silently substituted for the other.
- Fit each method's right-singular basis on its own calibration cycles under identical rules; at least eight calibration cycles per stratum; rank 3 is primary, ranks 1 and 5 are descriptive sensitivity checks. Do not subtract a learned mean waveform or add a learned template. Remove each row's scalar DC for projection and restore it for output. Undersupported groups bypass SVD equally for all methods and are flagged.
- TSA is arithmetic averaging of the test rows in each fixed stratum/channel. The half-cycle split uses chronological alternating rows, with the same pairs and equal counts for all arms. These are dependent within-record comparisons, not independent experiments.

## Evidence and failure checks

1. Show fixed S1/angle-bin-0 TSA waveforms and envelope TSA for the first test block. All methods share axes and acquisition-amplitude scale; no per-curve normalization.
2. Show the first eight Fourier-series amplitudes of the full-cycle TSA envelope. Frequency labels are h times the measured mean fp, NOT a 0.1 Hz continuous FFT. Averaging onto a nominated period guarantees a harmonic grid; its existence alone is not evidence of fault identity. Do not calculate a background SNR from the forced line grid or concatenate event windows to obtain a fault spectrum.
3. Report TSA output AC energy / mean input-row AC energy, and cross-validated prediction of original, unprojected held-out rows by the other half's processed TSA. The latter is a descriptive energy-explanation fraction, may be negative, and does not treat the projected test rows as ground truth.
4. Report within-condition harmonic-amplitude changes vs the unaligned arm with identical downstream processing, including all three blocks. BL checks whether apparently beneficial changes also occur without the target fault; no cross-condition physical amplitude comparison is justified by unverified sensitivities.
5. Validate the code on known phase-jittered pulse cycles, alternating-sign ringing, a known Fourier coefficient, frozen projection linearity, and held-out noninterference. Preserve all original data and prior results.

Scope: exploratory single previously inspected recording per condition; fixed signal band and inherited 21/31/84 teeth (geometry unresolved). This audit tests these TSA implementations, not every possible TSA/TPSVD choice or fault-contact ground truth.
