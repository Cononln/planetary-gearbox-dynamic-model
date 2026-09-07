function dx = pg_rhs_baseline(t, x, p)
%PG_RHS_BASELINE First-order state equation for the Stage-1 model.
% State x = [q; qdot], with carrier-angle-periodic M,C,K matrices.

n = p.map.n;
q = x(1:n);
qd = x(n+1:2*n);
phi_c = p.model.phi_c0 + p.kin.omega_c*t;
[M,C,K,meta] = pg_assemble_baseline(phi_c,p);
force = pg_baseline_force(t,phi_c,p,meta);
qdd = M\(force - C*qd - K*q);
dx = [qd; qdd];
end
