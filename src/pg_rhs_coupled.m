function dx = pg_rhs_coupled(t,x,p,fault,model)
%PG_RHS_COUPLED First-order state equation for the coupled 24-DOF model.

n = p.map.n+numel(model.frequencyHz);
q = x(1:n);
qd = x(n+1:2*n);
[M,C,K,meta] = pg_assemble_coupled(t,p,fault,model);
force = pg_coupled_force(t,p,meta);
qdd = M\(force-C*qd-K*q);
dx = [qd;qdd];
end
