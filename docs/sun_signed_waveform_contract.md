# Sun-fault signed-waveform trial

Use the available 50% sun-tooth-fault, nominal 600 rpm, L00 three-channel record.
This MAT file has chanvals (31,813,222 x 5 double), no sampling-rate field. Use the
51,200 Hz value from the existing MATLAB loader as an explicit inherited assumption;
columns 1-3 are vibration and 4/5 are encoder A/B. Encoder timing is a posterior
sanity check, not the source of the primary angular coordinate. No sensitivity
conversion or cross-record physical-amplitude comparison is performed.

Read source seconds 30-90. Estimate vibration motion as in the preceding pilot.
Use internal seconds 3-59, retaining complete sun-fault tidal rows only. For sun21,
ring84, three planets, the recurrence of faulty sun tooth + same planet + carrier
position is one carrier revolution (five sun rotations), with 12 aggregate sun
fault events and 84 mesh periods. This is not a repeat of the identity of every
healthy planet tooth. The inherited gear geometry remains provisional.

Each row uses 16,384 uniformly spaced angular samples. Compare no channel lag,
one whole matched-span signed Pearson lag, and one lag per sun tidal row. S1 is
fixed, S2/S3 are fitted independently, +/-410 original samples, no polarity flip,
no individual-event shifts, no extra within-channel inter-row registration.

Retain signed 4-10 kHz vibration throughout matching, SVD, row TSA and channel
fusion. The real FIR is mathematically the real part of the former complex
bandpass filter; no temporal envelope is computed. Original envelope-dependent
PSI is not used in this pilot. Fixed ranks 0 (bypass), 1, 2 are sensitivity controls.

Diagnostics: signed fusion RMS, channel coherent-energy fraction, descriptive
odd/even row-average correlation, direct Fourier spectra, and exploratory
sun-fault-spaced comb/background ratio. Use per-row spectral bins: carrier=1,
fault=12, mesh=84. State the much coarser ~2 Hz Fourier grid of a ~0.5 s TSA output;
averaging more rows improves repeatability but does not lengthen the output period.
The comb metric alone is not independent proof of fault specificity.

Figure contract: existing Python quantitative-grid style, 183 mm width, fixed
method colours, common amplitude scales, editable SVG/PDF and PNG/TIFF. Show one
complete sun tidal mean, not a repeated artificial 1 s trace; local zoom is chosen
from unchanged S1. Direct spectra and lag/row-consistency diagnostics carry the
remaining evidence. Export all plotted points and full arrays. No invented error
bars or significance: rows belong to one historical recording.
