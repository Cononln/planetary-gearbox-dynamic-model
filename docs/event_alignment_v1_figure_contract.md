# Event-alignment pilot figure contract

Backend: Python, retained from the user's established preference.
Archetype: quantitative grids with one event-local waveform panel as primary visual evidence.
Purpose: show whether event-constrained delays generalize to held-out sensors and improve signed pulse averaging relative to local waveform matching; negative findings remain visible.
Exports: 183 mm wide; Arial-compatible text; editable SVG/PDF, PNG preview and 600 dpi TIFF.
Panel plan: event evidence and candidate timing; accepted-event delay maps; common-mask raw/local/constrained pulse overlays; frozen-SVD TSA; whole-block envelope spectra; coverage and cross-block phase scores.
Statistics: one recording per condition, three disjoint time blocks; no independent-replicate p values or confidence intervals. Every plotted point must identify its condition, channel, event/angle stratum or time block in CSV/NPZ source data.
Integrity: examples chosen by fixed time and angle stratum, never highest measured gain. Same axes and event sets for before/after. Missing events remain explicit. Any y-axis scaling is shared within a comparison; no per-method amplitude scaling of measured output.
Reviewer risks: same-window adaptive matching; shared mesh interference; insufficient event support; signed phase differing after peak alignment; sparse full tidal repeats; previously inspected records; phase-only changes failing to improve fault/background contrast.

Post-run diagnostic figure: envelope-vs-signed coherence on exactly the frozen event shifts and common masks. Conclusion: constrained coarse shifts improve PF50 envelope concentration but do not establish signed waveform phase alignment. Two panels compare BL/PF50 using matched axes; three block-pair means are dependent summaries, not independent repetitions. Continuous-record Hilbert envelopes are calculated before event cropping. No new model tuning or fault-specificity claim is permitted from this diagnostic alone.
