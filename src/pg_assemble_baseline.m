function [M, C, K, meta] = pg_assemble_baseline(phi_c, p, meshStiffness)
%PG_ASSEMBLE_BASELINE Assemble M,C,K at a prescribed mean carrier angle.
% Stage 1 uses constant mean mesh stiffness; orientation is time periodic.

if nargin < 2
    p = pg_parameters;
end
if nargin < 3
    meshStiffness.sunPlanet = ...
        p.mesh.sunPlanet.kMean*ones(1,p.model.nPlanet);
    meshStiffness.ringPlanet = ...
        p.mesh.ringPlanet.kMean*ones(1,p.model.nPlanet);
end

M = pg_mass_matrix(p);
C = zeros(p.map.n);
K = zeros(p.map.n);

[C, K] = addGroundSupport(C, K, p.map.sun, p.support.sun);
[C, K] = addGroundSupport(C, K, p.map.ring, p.support.ring);
[C, K] = addGroundSupport(C, K, p.map.carrier, p.support.carrier);

% Planet-pin bearings couple planet translation to carrier translation/rotation.
pinVectors = cell(p.model.nPlanet, 2);
for i = 1:p.model.nPlanet
    psi = phi_c + p.gear.planetPhase0(i);
    bx = zeros(p.map.n,1);
    by = zeros(p.map.n,1);
    bx(p.map.planet{i}(1)) = 1;
    bx(p.map.carrier(1)) = -1;
    bx(p.map.carrier(3)) = p.gear.carrierRadius*sin(psi);
    by(p.map.planet{i}(2)) = 1;
    by(p.map.carrier(2)) = -1;
    by(p.map.carrier(3)) = -p.gear.carrierRadius*cos(psi);
    K = K + p.support.planetPin.kxy*(bx*bx.' + by*by.');
    C = C + p.support.planetPin.cxy*(bx*bx.' + by*by.');
    pinVectors{i,1} = bx;
    pinVectors{i,2} = by;
end

mesh = pg_mesh_vectors(phi_c, p);
meshDamping = repmat(struct('cSunPlanet',0,'cRingPlanet',0), ...
                          1, p.model.nPlanet);
Minv = diag(1./diag(M));
for i = 1:p.model.nPlanet
    bsp = mesh(i).bSunPlanet;
    brp = mesh(i).bRingPlanet;
    muSp = 1/(bsp.'*Minv*bsp);
    muRp = 1/(brp.'*Minv*brp);
    kSp = meshStiffness.sunPlanet(i);
    kRp = meshStiffness.ringPlanet(i);
    cSp = 2*p.mesh.sunPlanet.zeta*sqrt(max(kSp,0)*muSp);
    cRp = 2*p.mesh.ringPlanet.zeta*sqrt(max(kRp,0)*muRp);
    K = K + kSp*(bsp*bsp.') + kRp*(brp*brp.');
    C = C + cSp*(bsp*bsp.') + cRp*(brp*brp.');
    meshDamping(i).cSunPlanet = cSp;
    meshDamping(i).cRingPlanet = cRp;
end

% Remove numerical asymmetry introduced by floating-point accumulation.
K = (K + K.')/2;
C = (C + C.')/2;

meta.phi_c = phi_c;
meta.mesh = mesh;
meta.pinVectors = pinVectors;
meta.meshDamping = meshDamping;
meta.meshStiffness = meshStiffness;
end

function [C, K] = addGroundSupport(C, K, idx, support)
K(idx(1),idx(1)) = K(idx(1),idx(1)) + support.kxy;
K(idx(2),idx(2)) = K(idx(2),idx(2)) + support.kxy;
K(idx(3),idx(3)) = K(idx(3),idx(3)) + support.ktheta;
C(idx(1),idx(1)) = C(idx(1),idx(1)) + support.cxy;
C(idx(2),idx(2)) = C(idx(2),idx(2)) + support.cxy;
C(idx(3),idx(3)) = C(idx(3),idx(3)) + support.ctheta;
end
