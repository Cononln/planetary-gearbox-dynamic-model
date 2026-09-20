function [force,detail] = pg_pengyue_excitation_force(t,p,meta)
%PG_PENGYUE_EXCITATION_FORCE Sinusoidal TE plus balanced external torques.
% Eq. (2-28) uses e_v(t)=e_av*sin(omega_m*t+beta_v).

motion = pg_motion_state(t,p);
omegaMesh = motion.omega_mesh;
meshPhase = motion.meshPhaseRad;
force = zeros(p.map.n,1);
detail.sunPlanet = zeros(1,p.model.nPlanet);
detail.ringPlanet = zeros(1,p.model.nPlanet);

for ip = 1:p.model.nPlanet
    phSp = p.mesh.sunPlanet.tePhase(ip);
    phRp = p.mesh.ringPlanet.tePhase(ip);
    eSp = p.mesh.sunPlanet.teAmplitude*sin(meshPhase+phSp);
    eRp = p.mesh.ringPlanet.teAmplitude*sin(meshPhase+phRp);
    edSp = omegaMesh*p.mesh.sunPlanet.teAmplitude*cos(meshPhase+phSp);
    edRp = omegaMesh*p.mesh.ringPlanet.teAmplitude*cos(meshPhase+phRp);
    fSp = meta.meshStiffness.sunPlanet(ip)*eSp + ...
        meta.meshDamping(ip).cSunPlanet*edSp;
    fRp = meta.meshStiffness.ringPlanet(ip)*eRp + ...
        meta.meshDamping(ip).cRingPlanet*edRp;
    force = force + meta.mesh(ip).bSunPlanet*fSp + ...
        meta.mesh(ip).bRingPlanet*fRp;
    detail.sunPlanet(ip) = fSp;
    detail.ringPlanet(ip) = fRp;
end

% External drive/load terms. Tin is inferred by ideal power balance because
% the thesis states Tout=50 N m but does not report Tin for Fig. 3-14.
force(p.map.sun(3)) = force(p.map.sun(3))+p.operating.inputTorque;
force(p.map.carrier(3)) = force(p.map.carrier(3))-p.operating.outputTorque;
detail.inputTorque = p.operating.inputTorque;
detail.outputTorque = p.operating.outputTorque;
end
