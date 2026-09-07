# Fault-event arrival alignment V2: execution and figure contract

Objective: evaluate alignment of corresponding transient arrivals, not the signs of high-frequency ringing.
Backend: Python, using the established project runtime.
Data: reuse the previously inspected BL/PF50 records, existing vibration-derived motion and 4–10 kHz band. Development/calibration [3,15), validation [17,25); exploratory test blocks [27,37), [38,48), [49,59). No fresh-blind claim.
Algorithm: calibration onset-anchored envelope templates; sub-sample local envelope cross-correlation with multiple lag peaks; frozen relative periodic path; dynamic-programming common-event sequence; constraints select candidates rather than replace every measured lag with a smoothed prediction.
Onset proxy: 20% rising-edge crossing of a prominent envelope pulse, with local background correction and a bounded preceding-edge search. It is an operational sensor-response marker, not mechanical contact truth.
Primary evaluation: other-two-channel prediction of a fixed held-out-channel onset marker, on common supported events and with coverage. Compare a fixed-offset/median baseline and periodic joint inference. Known-arrival waveform controls verify timing accuracy separately.
Additional evaluation: residual inter-channel onset spread in a separate 0.8–2.5 kHz response band; envelope consistency and original-time-axis spectra; healthy and wrong-order controls. Cross-band markers can differ through dispersion and are corroborating proxies only.
Baselines: no correction; preceding waveform XCorr; local envelope XCorr; V1 constrained shifts; V2 joint envelope candidate selection. Main V2 comparison shares templates, event schedule, search limits and common evaluation masks with local envelope XCorr.
Selection: fixed small grid of path-constraint strengths, selected from validation other-channel prediction plus missingness penalty before reading V2 test metrics. Prior band selection is inherited and acknowledged.
Integrity: one record/state; three dependent block-pair summaries are not independent repetitions. Main plots before SVD, no inverse-amplitude compensation. Wide event displays retain the original pulse before shifting. No forced best-looking example.
Figures: quantitative comparison grid; fixed event packet overlays and means; original-axis envelope spectrum with fault-family/background values; healthy and wrong-order controls. 183 mm width, editable SVG, PNG preview, PDF and 600 dpi TIFF. Same scales and explicit sample/support definitions.
Decision: report improvements or failures as observed. Envelope consistency alone cannot establish correct fault identity. Any claimed speed/path separation remains limited by the common-phase gauge ambiguity.
