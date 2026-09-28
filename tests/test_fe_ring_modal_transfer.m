function report = test_fe_ring_modal_transfer
%TEST_FE_RING_MODAL_TRANSFER Asset-free checks of the FE path kernel.
% This synthetic fixture tests dimensions, periodic contact interpolation,
% zero-input response and agreement of the modal FRF with its equation.

model.frequencyHz = 1000;
model.omega = 2*pi*model.frequencyHz;
model.zeta = 0.02;
model.contactAngleRad = [0;pi];
model.inputShapeNormal = [1;1];
model.sensorShapeRadial = [1;2];
model.sensorAnglesLocalDeg = [0 90];
model.sensorAnglesGlobalDeg = [0 270];
model.sensorObservationStatus = {'synthetic 0';'synthetic 90'};

p.sensor.locationAngles = [0 pi/2];
p.gear.planetPhase0 = [0 2*pi/3 4*pi/3];
thetaFe = pg_fe_ring_contact_angle(0,p,model,1:3);
assert(max(abs(thetaFe-[0 4*pi/3 2*pi/3]))<1e-12, ...
    'Clockwise local-to-FE contact registration is inconsistent.');

psi = pg_fe_ring_input_shape([0;2*pi],model);
assert(isequal(size(psi),[2 1]) && max(abs(psi-1))<1e-12, ...
    'Contact interpolation is not periodic.');

f = 500;
[H,E] = pg_fe_ring_modal_path_gain(0,f,model);
w = 2*pi*f;
expected = -w^2/(model.omega^2-w^2+1i*2*model.zeta*model.omega*w);
assert(isequal(size(H),[1 2]) && abs(H(1,1)-expected)<1e-12, ...
    'FE modal FRF differs from the modal equation.');
assert(abs(H(1,2)-2*expected)<1e-12 && ...
    max(abs(E(:)-abs(H(:)).^2))<1e-12, ...
    'Sensor projection or energy is inconsistent.');

p.model.nPlanet = 3;
time = (0:255)'/16000;
force = zeros(numel(time),3);
carrier = zeros(numel(time),1);
out = pg_apply_fe_ring_transfer_causal(time,force,carrier,p,model);
assert(isequal(size(out.acceleration),[numel(time) 2]) && ...
    ~any(out.acceleration(:)), ...
    'Zero PR force should give zero modal sensor response.');

report = struct('passed',true,'modalDOF',1,'sensorCount',2, ...
    'sampleCount',numel(time));
fprintf('FE ring modal transfer smoke test PASSED.\n');
end
