# Signed-waveform pilot

Purpose: test the user's requested signed-waveform route with only the matching
representation and downstream fusion/FFT changed from the previous scope study.
No outcome is assumed. The study cannot establish fault-selective suppression.

- BL and PF50: one historical recording each; 4-10 kHz band waveforms at 51,200 Hz.
- Reference S1 is fixed. Compare unshifted, one whole-span signed Pearson lag, and
  one signed Pearson lag per full tidal row. S2/S3 are shifted relative to S1.
- Signed Pearson maximization over +/-410 samples; fractional interpolation;
  no polarity changes, circular wrap, individual-event shifts or extra row registration.
- Use the same vibration-derived motion coordinate, three 31-carrier-turn rows,
  and 524,288 points per row. The inherited 21/31/84 geometry is provisional.
- Independently recompute alignment, resampling, per-channel SVD and row TSA from
  stored band waveforms; compare new means with previously saved signed branches.
- Bypass SVD, rank 1, rank 2 are fixed sensitivity controls. Original PSI requires
  a modal envelope and is not used in this explicitly authorized pilot.
- Fuse three signed channel TSA waveforms arithmetically, then compute their direct
  Fourier series. No temporal rectification or analytic-envelope calculation.
  Magnitude of a Fourier coefficient is still a valid amplitude spectrum.
- This is angle-domain processing: frequency is equivalent frequency, k/mean row
  duration, not a direct uniformly time-sampled FFT. No claim of a 5.38 Hz
  low-frequency peak from a 4-10 kHz input is made.

Diagnostics: fused RMS, channel coherent-energy fraction, spectral mesh-line energy
fraction, and an exploratory fault-spaced comb/background ratio. The latter uses
three bins around multiples of 84 Fourier bins in 4-10 kHz, excluding +/-3 carrier
orders plus a margin around mesh lines. Matched background uses offsets +/-5..10
bins; both target and background use RMS amplitudes. This cannot distinguish a
fault from every other commensurate periodic component. Metrics are in-sample.

Figure contract: Python quantitative grids, 183 mm wide, Arial, editable SVG/PDF,
PNG/TIFF previews, shared scales and fixed colours. Waveforms show one second and
a short zoom around the strongest unchanged S1 sample in that second; this is a
response example, not confirmed fault-contact timing. Display-only thinning keeps
actual extrema, with every plotted value exported. Full numerical arrays are saved.
Spectra show 0-12 kHz and a fixed baseline-derived high-frequency zoom. No error
bars or significance claims: three dependent tidal rows, not independent records.

Acceptance: signed arithmetic cancellation unit test, exact FFT amplitude, lag
sign, full-rank SVD identity, no analytic-envelope call, cached signed-branch
agreement, same spans and fixed S1, source hashes unchanged, finite diagnostics,
Parseval check, and visual inspection of every figure.
