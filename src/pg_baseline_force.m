function [force, detail] = pg_baseline_force(t, phi_c, p, meta)
%PG_BASELINE_FORCE Deterministic transmission-error forcing interface.
% Stage 2 augments this interface with fault-induced TVMS forcing.

if nargin < 4
    [~,~,~,meta] = pg_assemble_baseline(phi_c, p);
end
omegaMesh = 2*pi*p.kin.f_mesh;
force = zeros(p.map.n,1);
detail.sunPlanet = zeros(1,p.model.nPlanet);
detail.ringPlanet = zeros(1,p.model.nPlanet);

for i = 1:p.model.nPlanet
    phSp = p.mesh.sunPlanet.tePhase(i);
    phRp = p.mesh.ringPlanet.tePhase(i);
    eSp = p.mesh.sunPlanet.teAmplitude*cos(omegaMesh*t + phSp);
    eRp = p.mesh.ringPlanet.teAmplitude*cos(omegaMesh*t + phRp);
    edSp = -omegaMesh*p.mesh.sunPlanet.teAmplitude*sin(omegaMesh*t + phSp);
    edRp = -omegaMesh*p.mesh.ringPlanet.teAmplitude*sin(omegaMesh*t + phRp);
    fSp = meta.meshStiffness.sunPlanet(i)*eSp + ...
          meta.meshDamping(i).cSunPlanet*edSp;
    fRp = meta.meshStiffness.ringPlanet(i)*eRp + ...
          meta.meshDamping(i).cRingPlanet*edRp;
    force = force + meta.mesh(i).bSunPlanet*fSp + ...
                    meta.mesh(i).bRingPlanet*fRp;
    detail.sunPlanet(i) = fSp;
    detail.ringPlanet(i) = fRp;
end
end
