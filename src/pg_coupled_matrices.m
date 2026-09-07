function [M,C,K,meta] = pg_coupled_matrices(phi_c,p,model,meshStiffness)
%PG_COUPLED_MATRICES Assemble bidirectionally coupled 18+6 DOF matrices.
% Ring-planet mesh deflection contains both rigid-body and local elastic-ring
% displacement, so each ring mesh generates off-diagonal coupling blocks.

if nargin < 3 || isempty(model)
    model = pg_ring_modal_model(p);
end
if nargin < 4
    meshStiffness.sunPlanet = ...
        p.mesh.sunPlanet.kMean*ones(1,p.model.nPlanet);
    meshStiffness.ringPlanet = ...
        p.mesh.ringPlanet.kMean*ones(1,p.model.nPlanet);
end

% Assemble supports, pins, and sun-planet meshes in the rigid subsystem.
rigidStiffness = meshStiffness;
rigidStiffness.ringPlanet(:) = 0;
[Mg,Cg,Kg,metaRigid] = pg_assemble_baseline(phi_c,p,rigidStiffness);

nq = p.map.n;
nr = numel(model.frequencyHz);
M = blkdiag(Mg,model.M);
C = blkdiag(Cg,model.C);
K = blkdiag(Kg,model.K);
Minv = diag(1./diag(M));
PsiInput = pg_ring_mode_shape(phi_c+p.gear.planetPhase0,model);

coupling = repmat(struct('g',zeros(nq+nr,1),'cMesh',0, ...
    'psiInput',zeros(nr,1)),1,p.model.nPlanet);
for i = 1:p.model.nPlanet
    g = [metaRigid.mesh(i).bRingPlanet; ...
         model.radialMeshProjection*PsiInput(i,:).'];
    kRp = meshStiffness.ringPlanet(i);
    muRp = 1/(g.'*Minv*g);
    cRp = 2*p.mesh.ringPlanet.zeta*sqrt(max(kRp,0)*muRp);
    K = K+kRp*(g*g.');
    C = C+cRp*(g*g.');
    coupling(i).g = g;
    coupling(i).cMesh = cRp;
    coupling(i).psiInput = PsiInput(i,:).';
end

M = (M+M.')/2;
C = (C+C.')/2;
K = (K+K.')/2;
meta.phi_c = phi_c;
meta.rigid = metaRigid;
meta.coupling = coupling;
meta.meshStiffness = meshStiffness;
meta.ringModel = model;
end
