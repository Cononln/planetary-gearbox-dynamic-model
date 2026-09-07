# V2 fault-frequency audit

Purpose: determine whether frozen event alignment increases planetary-fault spectral components, rather than inferring frequency-domain success from sharper event averages.
Backend: Python. Archetype: quantitative grids; 183 mm wide; editable SVG/PDF, PNG preview and 600 dpi TIFF.
Inputs: unchanged V2 event shifts and original cached band responses; previously inspected BL/PF50 records, three 10 s blocks per state. No refitting or parameter selection.
Physical interpretation: inherited teeth 21/31/84 remain provisional geometrically. Carrier frequency is mesh/84; a faulty planet's same-side recurrence is mesh/31. Twice that frequency may be strong when both contacts contribute but is not guaranteed to be the fundamental. Carrier modulation peaks must not be labeled as fault peaks.
Panel plan: actual time-domain envelope FFT in Hz, 0–50 Hz; fixed zooms of the first four fault harmonics; exact-order complex envelope coefficients and multi-channel spectral phasor retention.
Methods: unaligned, local envelope XCorr, joint envelope candidate matching. Freeze all shifts and use identical signal lengths/window/amplitude scale. Spectrum of the time-domain mean of individual channel envelopes is distinguished from averaging spectral magnitudes.
Metrics: retain original V2 fault-family/background definitions; add individual harmonic amplitudes and exact-order coherent-fusion retention |mean(Cj)|/mean(|Cj|). These are descriptive diagnostics, not new proposed indices or physical SNRs.
Statistics: three blocks from one record per state; no independent-repeat significance. Show or export every block, and label any averages. No zero padding marketed as better resolution (10 s time FFT bin spacing 0.1 Hz).
Integrity: continuous original-time-axis reconstruction, no event-window concatenation; no SVD, amplitude inverse compensation, or new fitting. Wide overview includes strong non-fault peaks. Post-hoc analysis and measured-response rather than source-force interpretation remain explicit.
