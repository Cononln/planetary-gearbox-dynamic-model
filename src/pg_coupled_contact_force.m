function detail = pg_coupled_contact_force(t,q,qd,p,meta)
%PG_COUPLED_CONTACT_FORCE Recover signed mesh forces in the 18+modal model.

nq = p.map.n;
nr = numel(meta.ringModel.frequencyHz);
motion = pg_motion_state(t,p);
omegaMesh = motion.omega_mesh;
meshPhaseRad = motion.meshPhaseRad;
detail.sunPlanet = zeros(1,p.model.nPlanet);
detail.ringPlanet = zeros(1,p.model.nPlanet);

for ip = 1:p.model.nPlanet
    phSp = p.mesh.sunPlanet.tePhase(ip);
    phRp = p.mesh.ringPlanet.tePhase(ip);
    eSp = p.mesh.sunPlanet.teAmplitude*cos(meshPhaseRad+phSp);
    eRp = p.mesh.ringPlanet.teAmplitude*cos(meshPhaseRad+phRp);
    edSp = -omegaMesh*p.mesh.sunPlanet.teAmplitude*sin(meshPhaseRad+phSp);
    edRp = -omegaMesh*p.mesh.ringPlanet.teAmplitude*sin(meshPhaseRad+phRp);

    gSp = [meta.rigid.mesh(ip).bSunPlanet;zeros(nr,1)];
    gRp = meta.coupling(ip).g;
    kSp = meta.meshStiffness.sunPlanet(ip);
    kRp = meta.meshStiffness.ringPlanet(ip);
    cSp = meta.rigid.meshDamping(ip).cSunPlanet;
    cRp = meta.coupling(ip).cMesh;

    detail.sunPlanet(ip) = ...
        kSp*(eSp-gSp.'*q)+cSp*(edSp-gSp.'*qd);
    detail.ringPlanet(ip) = ...
        kRp*(eRp-gRp.'*q)+cRp*(edRp-gRp.'*qd);
end
end
