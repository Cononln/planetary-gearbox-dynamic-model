function D = pg_ring_dynamic_stiffness(omega,model)
%PG_RING_DYNAMIC_STIFFNESS Complex modal dynamic-stiffness matrix.

D = model.K - omega^2*model.M + 1i*omega*model.C;
end
