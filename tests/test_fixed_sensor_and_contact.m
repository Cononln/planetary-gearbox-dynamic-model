function tests = test_fixed_sensor_and_contact
tests = functiontests(localfunctions);
end

function testDefaultSensorObservationIsBackwardCompatible(testCase)
p = pg_parameters;
model = pg_ring_modal_model(p);
C = pg_sensor_observation(p,model);
verifySize(testCase,C,[2,p.map.n+numel(model.frequencyHz)]);
verifyEqual(testCase,C(1,p.map.ring(1:2)),[1,0],'AbsTol',1e-12);
verifyEqual(testCase,C(2,p.map.ring(1:2)),[0,1],'AbsTol',1e-12);
end

function testSensorLocationAndSensitiveDirectionAreIndependent(testCase)
p = pg_parameters;
p.sensor.locationAngles = [0,2*pi/3,4*pi/3];
p.sensor.sensitiveAngles = [0,pi/2,pi];
model = pg_ring_modal_model(p);
C = pg_sensor_observation(p,model);
verifySize(testCase,C,[3,p.map.n+numel(model.frequencyHz)]);
verifyEqual(testCase,C(:,p.map.ring(1:2)), ...
    [1,0;0,1;-1,0],'AbsTol',1e-12);
end

function testFixedPathGeometryRepeatsEveryCarrierTurn(testCase)
p = pg_parameters;
g0 = pg_fixed_path_kinematics(0.37,p);
g1 = pg_fixed_path_kinematics(0.37+2*pi,p);
verifyEqual(testCase,g1.radialProjection,g0.radialProjection, ...
    'AbsTol',1e-12);
verifyEqual(testCase,g1.tangentialProjection,g0.tangentialProjection, ...
    'AbsTol',1e-12);
end

function testThreeSensorPathFrf(testCase)
p = pg_parameters;
p.sensor.locationAngles = [0,2*pi/3,4*pi/3];
p.sensor.sensitiveAngles = p.sensor.locationAngles;
model = pg_ring_modal_model(p);
path = pg_path_frf(p.kin.f_mesh,0.2,p,model);
verifySize(testCase,path.H,[3,p.model.nPlanet]);
verifyEqual(testCase,path.sensorAngles,p.sensor.locationAngles, ...
    'AbsTol',1e-12);
verifyTrue(testCase,all(isfinite(path.H),'all'));
end

function testContactLaws(testCase)
delta = [-2,-0.5,0.5,2];
deltaDot = zeros(size(delta));
[fLinear,~,~,closedLinear] = pg_contact_law( ...
    delta,deltaDot,10,0,'linear',1,true);
[fUni,~,~,closedUni] = pg_contact_law( ...
    delta,deltaDot,10,0,'unilateral',1,true);
[fBacklash,~,~,closedBacklash] = pg_contact_law( ...
    delta,deltaDot,10,0,'backlash',1,true);
verifyEqual(testCase,fLinear,10*delta,'AbsTol',1e-12);
verifyTrue(testCase,all(closedLinear));
verifyEqual(testCase,fUni,[0,0,0,10],'AbsTol',1e-12);
verifyEqual(testCase,closedUni,[false,false,false,true]);
verifyEqual(testCase,fBacklash,[-10,0,0,10],'AbsTol',1e-12);
verifyEqual(testCase,closedBacklash,[true,false,false,true]);
end

function testLinearContactStateMatchesAssembledEquation(testCase)
p = pg_parameters_loaded_v2;
p.contact.model = 'linear';
model = pg_ring_modal_model(p);
fault = pg_fault_case('none',0,p);
t = 0.013;
n = p.map.n+numel(model.frequencyHz);
q = linspace(-2e-7,2e-7,n).';
qd = linspace(3e-5,-3e-5,n).';
[~,C,K,meta] = pg_assemble_coupled(t,p,fault,model);
[force,~] = pg_coupled_force_loaded_v2(t,p,meta);
assembledResidual = C*qd+K*q-force;
state = pg_coupled_contact_state(t,q,qd,p,fault,model);
contactResidual = state.C0*qd+state.K0*q- ...
    state.operatingForce-state.contactForce;
verifyEqual(testCase,contactResidual,assembledResidual, ...
    'RelTol',1e-10,'AbsTol',1e-7);
end
