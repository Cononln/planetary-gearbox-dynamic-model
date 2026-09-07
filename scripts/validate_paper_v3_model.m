function report = validate_paper_v3_model
%VALIDATE_PAPER_V3_MODEL Fast mechanism checks for the paper-facing model.

setup_paths;
p = pg_parameters_paper_v3;
model = pg_ring_modal_model(p);
sensor = pg_sensor_configuration(p);
assert(sensor.nSensor==3);
assert(max(abs(rad2deg(sensor.locationAngles)-[0,120,240]))<1e-12);
assert(max(abs(sensor.locationAngles-sensor.sensitiveAngles))<1e-12);

% A controlled 1.5% carrier-speed fluctuation must move all coordinates
% through one common integrated phase trajectory.
pVariation = p;
pVariation.operating.speedProfile.enabled = true;
pVariation.operating.speedProfile.fraction = 0.015;
pVariation.operating.speedProfile.frequencyHz = 0.65;
time = (0:1/51200:0.8)';
motion = pg_motion_state(time,pVariation);
numericalSpeed = gradient(motion.phi_c,time);
speedRmse = sqrt(mean((numericalSpeed-motion.omega_c).^2));
assert(speedRmse<2e-5);

fault = pg_fault_case('sun',0.25,pVariation);
[~,detail] = pg_tvms(time,pVariation,fault);
loss = detail.faultStiffnessLossSunPlanet+ ...
    detail.faultStiffnessLossRingPlanet;
truth = pg_fault_event_truth(time,loss,pVariation,fault);
assert(truth.completeEventCount>=15);
eventFrequencyHz = 1/mean(diff(truth.eventTime));
assert(abs(eventFrequencyHz-p.kin.f_sun_fault)<1.0);
assert(all(truth.peakStiffnessLossNPerM>0));

frequencyHz = (1700:25:2050)';
carrierAngle = linspace(0,2*pi,73);
pathTruth = pg_multichannel_phase_truth_field( ...
    frequencyHz,carrierAngle,pVariation,model,1,1);
assert(isequal(size(pathTruth.H), ...
    [numel(frequencyHz),numel(carrierAngle),3]));
assert(all(isfinite(pathTruth.H(:))));
assert(max(abs(pathTruth.H(:,1,:)-pathTruth.H(:,end,:)),[],'all')<1e-9);

report.sensorAnglesDeg = sensor.locationAnglesDeg;
report.speedProfileFraction = pVariation.operating.speedProfile.fraction;
report.speedProfileFrequencyHz = pVariation.operating.speedProfile.frequencyHz;
report.phaseSpeedRmseRadS = speedRmse;
report.completeFaultEventCount = truth.completeEventCount;
report.faultEventFrequencyHz = eventFrequencyHz;
report.nominalFaultFrequencyHz = p.kin.f_sun_fault;
report.pathTruthPeriodicityResidual = max(abs(pathTruth.H(:,1,:)- ...
    pathTruth.H(:,end,:)),[],'all');
report.pathTruthSensorCount = size(pathTruth.H,3);

summary = struct2table(report);
writetable(summary,fullfile('results','paper_v3_model_validation.csv'));
fprintf('Paper-v3 model validation PASSED.\n');
fprintf('  sensors (deg)             : %.0f / %.0f / %.0f\n', ...
    report.sensorAnglesDeg);
fprintf('  speed trajectory RMSE     : %.3e rad/s\n', ...
    report.phaseSpeedRmseRadS);
fprintf('  event frequency / nominal : %.3f / %.3f Hz\n', ...
    report.faultEventFrequencyHz,report.nominalFaultFrequencyHz);
end
