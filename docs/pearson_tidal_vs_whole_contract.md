# Pearson lag search: full tidal segments versus whole record

Question: does a separate channel lag per complete tidal repeat improve alignment and downstream averaging over a single lag for the same complete analyzed record?

Interpretation: within each tidal period, align S2 and S3 to S1; compare with one S2/S3-to-S1 lag over the concatenated analyzed time span. This is not same-channel alignment of successive tidal rows. Reference S1 is fixed, not selected for the best outcome.

Pearson plus cross-correlation means lag-dependent, locally mean/variance-normalized cross-correlation. Maximize signed Pearson, not absolute correlation. Primary search uses the 4–10 kHz envelope with the already established 0.12 ms smoothing. Signed band-waveform search is a separately labelled comparison. Search +/-8 ms in both scopes, with the same preprocessing and no event-specific adjustment.

Geometry: inherited ring84/planet31 yields complete recurrence after 31 carrier turns = 84 same-side fault cycles. Begin at 3 s, derive exact period boundaries from the existing vibration-estimated phase, and retain all three full tidal repeats before 59 s. Whole-record fitting uses exactly the same complete span; partial tail is excluded for both. Parameter geometry is provisional, as in previous analyses.

Downstream: build genuine full-tidal rows on a shared angle grid, 524,288 points per row; confirm equivalent angular sample rate stays above 22 kHz for the 4–10 kHz band. Compare no alignment, one whole-record lag, and one lag per tidal period. In each case use the same rows, SVD rank 1 or 2 (rank0 is TSA alone), and three-cycle TSA. Three-row rank3 would be an identity, not denoising. This new protocol deliberately differs from the previous frozen-calibration local/fault-cycle SVD and must not be pooled with it.

Representations: signed vibration SVD/TSA followed by channel-envelope averaging, and envelope SVD/TSA followed by channel averaging. No polarity flipping or inverse amplitude correction. No circular time-shift wraparound: read padded original samples using windowed-sinc interpolation. The remaining common vibration-derived angle resampling is identical across arms.

Evidence: plot the measured S2/S3 lags for the three full periods, within-period Pearson before/after (explicitly in-sample), matched cross-band event-marker spread as a response proxy, actual full-tidal TSA waveforms/envelopes, full-tidal Fourier-series spectra and fp harmonics. Fourier grid spacing ~1/15.6 Hz is an equivalent angle-domain frequency scale, not a zero-padded or direct uniformly sampled time FFT. Record all source values, per-period counts, rank sensitivity, BL control and search-boundary hits.

Risks: 15.6 s is itself much longer than carrier path variation. One lag per tidal period does not compensate path variation within that period. Optimized Pearson is not proof of aligning fault impulses. Only three complete rows and one historical record per state are available. SVD and alignment are in-sample descriptive processing here, not a blind holdout or independent significance test. The same time support, rank and fusion prevent attributing downstream differences to inconsistent averaging.

Figures: Python only, nature-figure quantitative grids, shared scales and stable method colours; 183 mm widths, editable SVG/PDF and PNG/TIFF previews. Panel conclusions will follow the computed results. Export numerical data, numerical checks and visual QA with the report.
