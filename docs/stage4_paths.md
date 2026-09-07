# Stage 4 — Coupled ring response and fixed dual-sensor paths

The 18 rigid-body coordinates and six analytical ring coordinates form a
24-DOF second-order system. For ring-planet mesh i, the elastic approach is

delta_rp_i = b_rp_i^T q_g
             + [-sin(alpha_rp)] psi(theta_i)^T eta
             - e_rp_i.

The resulting 24-entry projection vector generates k*g_i*g_i^T and
c*g_i*g_i^T. Its off-diagonal blocks couple the gear and ring subsystems in
both directions. This is stronger than applying a post-processing amplitude
window to an already generated signal.

The two observation rows are fixed in the housing frame:

- sensor 0 is at theta=0 and measures radial outward acceleration (+x);
- sensor 90 is at theta=pi/2 and measures radial outward acceleration (+y).

Each observation includes both rigid ring translation and local elastic modal
acceleration. Both sensors share the same acquisition clock; no electronic
clock offset is added.

For isolating the physical path, the harmonic path FRF uses healthy mean mesh
stiffness and a unit equivalent mesh-approach force at one moving ring mesh:

H_j_i(omega,phi_c) = -omega^2 C_j
                     [-omega^2 M + i omega C + K]^-1 g_i.

The three source paths are returned separately. Fault-force amplitudes can
later weight them without redefining the propagation path itself.
