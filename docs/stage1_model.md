# Stage 1 — 18-DOF lumped-parameter model

The generalized coordinate is

q = [q_s, q_r, q_c, q_p1, q_p2, q_p3]^T,

where every body coordinate is q_g = [x_g, y_g, theta_g]^T.

The baseline equation is

M qdd + C(phi_c) qd + K(phi_c) q = F_TE(t, phi_c).

The carrier mean angle is phi_c = omega_c t. The matrices change because the
three line-of-action vectors and the three planet-pin directions rotate past
the fixed inertial x/y frame. The dynamic carrier theta coordinate is retained
as a small torsional perturbation.

For planet i at psi_i = phi_c + 2*pi*(i-1)/3, the external sun-planet normal is

n_sp = sin(alpha_sp)*e_r + cos(alpha_sp)*e_t,

and the internal ring-planet normal is

n_rp = -sin(alpha_rp)*e_r + cos(alpha_rp)*e_t.

The mesh deflections are delta_sp_i = b_sp_i^T q - e_sp_i and
delta_rp_i = b_rp_i^T q - e_rp_i. Each mesh adds k*b*b^T and c*b*b^T to the
system matrices. Mesh damping is generated from a specified damping ratio and
the projection-vector effective mass.

The planet-pin relative translations also contain carrier torsional motion:

delta_pin_x = x_pi - x_c + r_c*sin(psi_i)*theta_c,

delta_pin_y = y_pi - y_c - r_c*cos(psi_i)*theta_c.

The ring keeps all three coordinates but has a strong housing support, so its
mean rotation is fixed without deleting the ring response required by Stage 3.

Bearing stiffness, damping, and constant mean mesh stiffness are currently
calibration parameters. CAD-derived mass, inertia, geometry, tooth numbers,
working pressure angle, and operating kinematics are kept separate from those
estimated parameters in pg_parameters.m.
