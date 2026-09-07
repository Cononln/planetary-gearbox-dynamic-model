function [sensorForce,weights] = pg_lpm_path_force(t,meshForce,p)
%PG_LPM_PATH_FORCE Map moving ring-mesh sources to fixed sensor sites.
% A periodic von-Mises proximity window is used instead of an arbitrary
% finite visibility cut-off.  The nonzero floor preserves the stationary
% mesh-frequency component, while higher spatial harmonics create the
% physically expected carrier/planet-passing modulation sidebands.

phiCarrier = p.model.phi_c0+p.kin.omega_c*t;
planetAngle = phiCarrier+p.gear.planetPhase0;
sensor = pg_sensor_configuration(p);
nSensor = sensor.nSensor;
weights = zeros(nSensor,p.model.nPlanet);
sensorForce = zeros(nSensor,1);

for is = 1:nSensor
    delta = planetAngle-sensor.locationAngles(is);
    proximity = exp(p.sensor.pathKappa*(cos(delta)-1));
    weights(is,:) = p.sensor.pathFloor+ ...
        (1-p.sensor.pathFloor)*proximity;
    sourceForce = p.sensor.planetSourceGain(:).* ...
        meshForce.ringPlanet(:);
    sensorForce(is) = weights(is,:)*sourceForce;
end
end
