# Stage 5 — Phase-difference and group-delay truth fields

For a specified moving ring-planet mesh source, the two fixed-sensor
acceleration transfer functions are H0(f,phi_c) and H90(f,phi_c). The truth
fields are defined by

Delta phi_true(f,phi_c) = arg[H90(f,phi_c) H0*(f,phi_c)],

Delta tau_g_true(f,phi_c)
    = -partial Delta phi_unwrapped(f,phi_c) / partial omega.

Phase is unwrapped along frequency before differentiation. Near a transfer
zero, phase and its derivative are not identifiable in measured data.
Accordingly, the output includes a reliability field and a valid mask derived
from the two path magnitudes and their amplitude balance. The unmasked raw
field is still retained in the MAT file.

Three angle-only summaries are computed as path-strength-weighted group delay
inside the 1.83, 3.03, and 4.98 kHz modal bands. They are periodic functions of
carrier angle, not a single cross-correlation lag.

The validation creates the same known broadband excitation behind both paths.
It compares:

1. no phase correction;
2. the full frequency correction measured at carrier angle zero and held
   static for the complete revolution;
3. the complete angle-frequency truth correction.

The third correction must reduce weighted cross-phase residual to numerical
precision. The static reference correction remains nonzero because it cannot
represent the moving physical path.

Primary outputs:

- results/phase_truth_field.mat: complete complex FRFs, phase, group delay,
  reliability, valid mask, and band summaries;
- results/phase_truth_band_delay.csv: encoder-angle lookup table for three
  resonance bands;
- results/step5_phase_truth.png: manuscript-oriented diagnostic overview.
