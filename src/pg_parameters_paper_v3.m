function p = pg_parameters_paper_v3
%PG_PARAMETERS_PAPER_V3 Paper-facing three-channel loaded model settings.
% This keeps the legacy two-channel examples reproducible while providing
% one explicit configuration for the actual 0/120/240-degree experiment.

p = pg_parameters_loaded_v2;
p.model.name = 'ZLS160 three-channel loaded 18+7 DOF effective-path v3.0';

% The three accelerometers are mounted at the stated circumferential
% positions.  Their sensitive axes are taken as outward radial directions;
% change only sensitiveAngles if the experimental mounting directions are
% later measured to be oblique.
p.sensor.locationAngles = [0,2*pi/3,4*pi/3];
p.sensor.sensitiveAngles = p.sensor.locationAngles;
p.sensor.angles = p.sensor.locationAngles;
p.sensor.mountingDescription = [ ...
    'Three fixed radial accelerometers at 0, 120, and 240 degrees.'];

% Controlled speed variation for the tacholess/event-alignment benchmark.
% It is disabled for the nominal response.  When enabled, pg_motion_state
% applies the same trajectory to every moving-coordinate term.
p.operating.speedProfile.enabled = false;
p.operating.speedProfile.fraction = 0;
p.operating.speedProfile.frequencyHz = 0.65;
p.operating.speedProfile.phaseRad = 0;
p.operating.speedProfile.definition = [ ...
    'omega_c(t)=omega_c0[1+a cos(2 pi f_v t+phi_v)]'];

p.v3.assumption = [ ...
    'Three-channel effective transfer-path model. Modal frequencies are ' ...
    'used as a low-order path prior; absolute FRF amplitude and phase are ' ...
    'not claimed without measured mode shapes or hammer-FRF calibration.'];
p.v3.primaryFaultModel = [ ...
    'Partial face-width broken tooth removed to the root circle; ' ...
    'severity is missing face-width fraction.'];
end
