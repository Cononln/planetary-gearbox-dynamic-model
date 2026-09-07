function Cacc = pg_lpm_sensor_observation(p)
%PG_LPM_SENSOR_OBSERVATION Project rigid-ring acceleration to fixed sensors.
% The local path-observer acceleration is added separately by the simulator.

sensor = pg_sensor_configuration(p);
Cacc = zeros(sensor.nSensor,p.map.n);
for is = 1:sensor.nSensor
    theta = sensor.sensitiveAngles(is);
    er = [cos(theta),sin(theta)];
    Cacc(is,p.map.ring(1:2)) = er;
end
end
