function Cacc = pg_sensor_observation(p,model)
%PG_SENSOR_OBSERVATION Acceleration observation rows for fixed radial sensors.
% Rigid-ring translation is projected on each accelerometer sensitive axis.
% Ring modal coordinates are radial, so their contribution also includes
% the projection between the local radial direction and that axis.

nq = p.map.n;
nr = numel(model.frequencyHz);
sensor = pg_sensor_configuration(p);
Cacc = zeros(sensor.nSensor,nq+nr);
for j = 1:sensor.nSensor
    thetaLocation = sensor.locationAngles(j);
    thetaSensitive = sensor.sensitiveAngles(j);
    er = [cos(thetaSensitive),sin(thetaSensitive)];
    Cacc(j,p.map.ring(1:2)) = er;
    radialProjection = cos(thetaSensitive-thetaLocation);
    Cacc(j,nq+(1:nr)) = radialProjection* ...
        pg_ring_observation(thetaLocation,model);
end
end
