function dx = pg_rhs(t,x,p,fault)
%PG_RHS First-order state equation with time-varying mesh stiffness.

n = p.map.n;
q = x(1:n);
qd = x(n+1:2*n);
phi_c = p.model.phi_c0 + p.kin.omega_c*t;
[M,C,K,meta] = pg_assemble_system(t,p,fault);
force = pg_baseline_force(t,phi_c,p,meta);
qdd = M\(force-C*qd-K*q);
dx = [qd;qdd];
end
