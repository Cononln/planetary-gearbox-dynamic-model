function detail = pg_mesh_contact_force(t,q,qd,p,meta)
%PG_MESH_CONTACT_FORCE Recover physical mesh forces from the LPM solution.
% For delta=b'*q-e, the generalized equation is
% M*qdd+C*qdot+K*q = b*(k*e+c*edot).  Hence the signed contact force along
% each line of action is k*(e-b'*q)+c*(edot-b'*qdot).

omegaMesh = 2*pi*p.kin.f_mesh;
detail.sunPlanet = zeros(1,p.model.nPlanet);
detail.ringPlanet = zeros(1,p.model.nPlanet);

for ip = 1:p.model.nPlanet
    phSp = p.mesh.sunPlanet.tePhase(ip);
    phRp = p.mesh.ringPlanet.tePhase(ip);
    eSp = p.mesh.sunPlanet.teAmplitude*cos(omegaMesh*t+phSp);
    eRp = p.mesh.ringPlanet.teAmplitude*cos(omegaMesh*t+phRp);
    edSp = -omegaMesh*p.mesh.sunPlanet.teAmplitude*sin(omegaMesh*t+phSp);
    edRp = -omegaMesh*p.mesh.ringPlanet.teAmplitude*sin(omegaMesh*t+phRp);

    bSp = meta.mesh(ip).bSunPlanet;
    bRp = meta.mesh(ip).bRingPlanet;
    kSp = meta.meshStiffness.sunPlanet(ip);
    kRp = meta.meshStiffness.ringPlanet(ip);
    cSp = meta.meshDamping(ip).cSunPlanet;
    cRp = meta.meshDamping(ip).cRingPlanet;

    detail.sunPlanet(ip) = kSp*(eSp-bSp.'*q)+cSp*(edSp-bSp.'*qd);
    detail.ringPlanet(ip) = kRp*(eRp-bRp.'*q)+cRp*(edRp-bRp.'*qd);
end
end
