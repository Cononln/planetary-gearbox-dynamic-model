function [B,thetaMesh] = pg_ring_input_matrix(phi_c,p,model)
%PG_RING_INPUT_MATRIX Map three moving ring-mesh radial forces to modal force.

thetaMesh = phi_c + p.gear.planetPhase0;
PsiInput = pg_ring_mode_shape(thetaMesh,model);
B = model.radialMeshProjection*PsiInput.';
end
