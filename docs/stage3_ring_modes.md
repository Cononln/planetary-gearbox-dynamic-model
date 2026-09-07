# Stage 3 — Frequency-calibrated analytical ring modes

Only modal frequencies were supplied; no finite-element node displacement
vectors or original result database exists. The spatial basis is therefore an
analytical circular-ring surrogate and must not be described as an FE modal
reduction.

The three near-degenerate pairs below are used in the 0–5 kHz band:

| Supplied mode | Frequency (Hz) | Analytical radial shape |
|---|---:|---|
| 7 | 1832.565469 | cos(2 theta) |
| 8 | 1845.781068 | sin(2 theta) |
| 9 | 3010.794375 | cos(3 theta) |
| 10 | 3053.891626 | sin(3 theta) |
| 11 | 4976.103825 | cos(4 theta) |
| 12 | 4976.220685 | sin(4 theta) |

The actual supplied frequency of every member is retained rather than forcing
an exactly degenerate pair. This small splitting is important: the cosine and
sine members then acquire different dynamic phases near resonance, producing
a carrier-angle-dependent point-to-point phase response.

With radial modes normalized to unit peak displacement, the generalized mass
of each non-axisymmetric sine/cosine mode is approximated as half the CAD ring
mass. Modal damping is initially 1% and remains a calibration parameter.

Three ring-planet mesh locations rotate with the carrier. Their angles are
theta_i = phi_c + 2*pi*(i-1)/3. A radial point force at theta_i contributes
the generalized modal force psi_r(theta_i)*F_i. The ring-planet line-of-action
force is projected radially using -sin(alpha_rp).
