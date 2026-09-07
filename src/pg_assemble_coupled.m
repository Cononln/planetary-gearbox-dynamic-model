function [M,C,K,meta] = pg_assemble_coupled(t,p,fault,model)
%PG_ASSEMBLE_COUPLED Time-domain 24-DOF assembly with Stage-2 TVMS.

if nargin < 4 || isempty(model)
    model = pg_ring_modal_model(p);
end
motion = pg_motion_state(t,p);
phi_c = motion.phi_c;
[meshStiffness,tvmsDetail] = pg_tvms(t,p,fault);
[M,C,K,meta] = pg_coupled_matrices(phi_c,p,model,meshStiffness);
meta.tvms = tvmsDetail;
meta.fault = fault;
end
