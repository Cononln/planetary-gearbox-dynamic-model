function path = pg_path_frf(frequencyHz,phi_c,p,model)
%PG_PATH_FRF Acceleration FRF from each moving ring mesh to fixed sensors.
% Healthy mean stiffness is used so the result isolates the transfer path
% from fault-force amplitude modulation.

meanStiffness.sunPlanet = ...
    p.mesh.sunPlanet.kMean*ones(1,p.model.nPlanet);
meanStiffness.ringPlanet = ...
    p.mesh.ringPlanet.kMean*ones(1,p.model.nPlanet);
[M,C,K,meta] = pg_coupled_matrices(phi_c,p,model,meanStiffness);
Cacc = pg_sensor_observation(p,model);
sensor = pg_sensor_configuration(p);
omega = 2*pi*frequencyHz;
D = K-omega^2*M+1i*omega*C;
path.H = complex(zeros(sensor.nSensor,p.model.nPlanet));
for i = 1:p.model.nPlanet
    path.H(:,i) = -omega^2*Cacc*(D\meta.coupling(i).g);
end
path.sensorAngles = sensor.locationAngles;
path.sensorSensitiveAngles = sensor.sensitiveAngles;
path.meshAngles = phi_c+p.gear.planetPhase0;
path.fixedGeometry = pg_fixed_path_kinematics(phi_c,p);
path.frequencyHz = frequencyHz;
path.phi_c = phi_c;
end
