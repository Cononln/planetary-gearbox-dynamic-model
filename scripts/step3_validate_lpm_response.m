function report = step3_validate_lpm_response
%STEP3_VALIDATE_LPM_RESPONSE Validate the default sensor-response chain.

setup_paths;
p = pg_parameters;
healthy = pg_fault_case('none',0,p);
severe = pg_fault_case('sun',0.50,p);
duration = 0.25;
discard = 0.10;

rHealthy = pg_simulate_lpm_response(duration,p,healthy);
rSevere = pg_simulate_lpm_response(duration,p,severe);
keep = rHealthy.time>=discard;
aHealthy = rHealthy.acceleration(keep,:);
aSevere = rSevere.acceleration(keep,:);
rmsHealthy = sqrt(mean(aHealthy.^2,1));
rmsSevere = sqrt(mean(aSevere.^2,1));

assert(p.map.n==18);
assert(isequal(size(rHealthy.acceleration,2),2));
assert(all(isfinite(rHealthy.acceleration),'all'));
assert(all(isfinite(rSevere.acceleration),'all'));
assert(all(rmsSevere>rmsHealthy));
assert(abs(mean(p.sensor.planetSourceGain)-1)<1e-12);
assert(range(p.sensor.planetSourceGain)>0);

dummy.sunPlanet = ones(1,p.model.nPlanet);
dummy.ringPlanet = ones(1,p.model.nPlanet);
[~,w0] = pg_lpm_path_force(0,dummy,p);
[~,wPeriod] = pg_lpm_path_force(1/p.kin.f_c,dummy,p);
periodicResidual = norm(wPeriod-w0,'fro')/norm(w0,'fro');
assert(periodicResidual<1e-12);
assert(norm(w0(1,:)-w0(2,:))>1e-3);

report.nGearDof = p.map.n;
report.sensorAnglesDeg = rad2deg(p.sensor.angles);
report.rmsHealthy = rmsHealthy;
report.rmsSunFault50 = rmsSevere;
report.pathPeriodicityResidual = periodicResidual;
report.pathNaturalFrequencyHz = p.sensor.pathNaturalFrequencyHz;
report.modelType = rHealthy.modelType;

fprintf('Stage 3 default LPM response PASSED\n');
fprintf('  gearbox DOFs                 : %d\n',report.nGearDof);
fprintf('  sensor angles                : %.0f / %.0f deg\n', ...
    report.sensorAnglesDeg);
fprintf('  healthy RMS 0/90             : %.4g / %.4g m/s^2\n', ...
    rmsHealthy);
fprintf('  50%% sun-fault RMS 0/90       : %.4g / %.4g m/s^2\n', ...
    rmsSevere);
fprintf('  path periodicity residual    : %.3e\n',periodicResidual);
end
