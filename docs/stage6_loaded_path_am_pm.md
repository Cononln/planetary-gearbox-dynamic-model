# Stage 6 — Loaded ring-path amplitude and phase modulation

The loaded 18+7-DOF model separates the mesh excitation from the moving
ring-to-sensor propagation path. For sensor `s` and planet mesh `p`, the
frozen-angle complex acceleration path is

    H_s,p(f,phi_c) = -(2*pi*f)^2 C_s
                     [K(phi_c)-omega^2 M+i omega C(phi_c)]^-1 g_p(phi_c).

The three paths are stored separately. For the healthy coherent case, their
in-phase sum is used. Other complex source weights can represent unequal
planet loading or a source localized at one planet without redefining the
structural path.

For a narrowband carrier at `f0`, the angle-periodic complex path can be
written as

    H_s(phi_c) = A_s(phi_c) exp(i psi_s(phi_c)).

`A_s` produces amplitude modulation and `psi_s` produces phase modulation.
Because `phi_c = 2*pi*f_c*t`, the Fourier-series coefficient of order `k`
appears at `f0 + k*f_c`. The analysis decomposes the full complex path into
AM-only and PM-only reconstructions, then compares their predicted orders
with the healthy time-domain simulation.

Two carriers are used:

- 168 Hz, the fundamental mesh center;
- 1848 Hz, the 11th mesh harmonic inside the first elastic-ring band.

The dual-channel path truth is

    Delta psi(f,phi_c) = arg[H_90(f,phi_c) H_0*(f,phi_c)].

This angle-frequency field is the model truth needed by a later path-aware
phase-alignment method. It must not be replaced by one constant time lag.

Primary outputs:

- `results/loaded_path_am_pm.png`;
- `results/loaded_single_planet_path_am_pm.png`;
- `results/loaded_path_sideband_decomposition.png`;
- `results/loaded_path_truth_maps.png`;
- `results/loaded_path_phase_lookup.csv` (1024 encoder angles, three paths);
- `results/loaded_path_am_pm_validation.csv` (path/sideband agreement and
  phase-lookup reliability);
- `results/loaded_path_am_pm_metrics.csv`;
- `results/loaded_path_sideband_decomposition.csv`;
- `results/loaded_path_am_pm.mat`.
