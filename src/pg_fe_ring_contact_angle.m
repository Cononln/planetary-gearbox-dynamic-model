function thetaGlobal = pg_fe_ring_contact_angle(carrierAngleRad,p,model,planetIndex)
%PG_FE_RING_CONTACT_ANGLE Register source-model contact angle to FE axes.
% Source local 0 deg is the SENSOR_0 direction. The exported sensor labels
% increase clockwise, while FE atan2 angles increase counterclockwise.
% The same transform is used for modal force input and event diagnostics.

thetaLocal = carrierAngleRad + p.gear.planetPhase0(planetIndex) - ...
    p.sensor.locationAngles(1);
thetaGlobal = mod(deg2rad(model.sensorAnglesGlobalDeg(1)) - thetaLocal,2*pi);
end
