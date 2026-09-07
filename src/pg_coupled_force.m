function [force,detail] = pg_coupled_force(t,p,meta)
%PG_COUPLED_FORCE Transmission-error forcing for the 24-DOF system.

nq = p.map.n;
nr = numel(meta.ringModel.frequencyHz);
force = zeros(nq+nr,1);
motion = pg_motion_state(t,p);
omegaMesh = motion.omega_mesh;
meshPhaseRad = motion.meshPhaseRad;
detail.sunPlanet = zeros(1,p.model.nPlanet);
detail.ringPlanet = zeros(1,p.model.nPlanet);

for i = 1:p.model.nPlanet
    phSp = p.mesh.sunPlanet.tePhase(i);
    phRp = p.mesh.ringPlanet.tePhase(i);
    eSp = p.mesh.sunPlanet.teAmplitude*cos(meshPhaseRad+phSp);
    eRp = p.mesh.ringPlanet.teAmplitude*cos(meshPhaseRad+phRp);
    edSp = -omegaMesh*p.mesh.sunPlanet.teAmplitude*sin(meshPhaseRad+phSp);
    edRp = -omegaMesh*p.mesh.ringPlanet.teAmplitude*sin(meshPhaseRad+phRp);

    kSp = meta.meshStiffness.sunPlanet(i);
    kRp = meta.meshStiffness.ringPlanet(i);
    cSp = meta.rigid.meshDamping(i).cSunPlanet;
    cRp = meta.coupling(i).cMesh;
    fSp = kSp*eSp+cSp*edSp;
    fRp = kRp*eRp+cRp*edRp;
    gSp = [meta.rigid.mesh(i).bSunPlanet;zeros(nr,1)];
    force = force+gSp*fSp+meta.coupling(i).g*fRp;
    detail.sunPlanet(i) = fSp;
    detail.ringPlanet(i) = fRp;
end
end
