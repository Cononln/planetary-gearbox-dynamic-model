function [M,C,K,meta] = pg_assemble_system(t,p,fault)
%PG_ASSEMBLE_SYSTEM Assemble the time-periodic system including Stage-2 TVMS.

phi_c = p.model.phi_c0 + p.kin.omega_c*t;
[stiffness,tvmsDetail] = pg_tvms(t,p,fault);
[M,C,K,meta] = pg_assemble_baseline(phi_c,p,stiffness);
meta.tvms = tvmsDetail;
meta.fault = fault;
end
