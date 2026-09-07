function geometry = pg_fixed_path_kinematics(phi_c,p)
%PG_FIXED_PATH_KINEMATICS Moving mesh geometry in the fixed sensor frame.
% This is the explicit coordinate layer used in the transmission-path
% model. The gearbox DOFs already live in an inertial frame, so the sensor
% observation matrix is fixed; carrier-periodic path variation enters via
% the moving mesh/source angles and the dynamic matrices.

if ~isscalar(phi_c) || ~isfinite(phi_c)
    error('pg_fixed_path_kinematics:CarrierAngle', ...
        'phi_c must be one finite carrier angle in radians.');
end
sensor = pg_sensor_configuration(p);
meshAngles = phi_c+p.gear.planetPhase0(:).';

geometry.phi_c = phi_c;
geometry.carrierToFixed = [cos(phi_c),-sin(phi_c); ...
                           sin(phi_c), cos(phi_c)];
geometry.meshAngles = meshAngles;
geometry.meshToSensorAngle = sensor.locationAngles(:)-meshAngles;
geometry.radialProjection = cos(sensor.sensitiveAngles(:)-meshAngles);
geometry.tangentialProjection = sin(sensor.sensitiveAngles(:)-meshAngles);
geometry.sensorLocationAngles = sensor.locationAngles;
geometry.sensorSensitiveAngles = sensor.sensitiveAngles;
end
