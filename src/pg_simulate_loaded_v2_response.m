function response = pg_simulate_loaded_v2_response(duration,p,fault,model,options)
%PG_SIMULATE_LOADED_V2_RESPONSE Loaded rigid-plus-modal ring simulation.
% Constant drive/load torques establish mesh preload.  A loaded equilibrium
% is used as the initial state, and the time-varying mesh stiffness then
% generates the dynamic response.

arguments
    duration (1,1) double {mustBePositive}
    p struct
    fault struct
    model struct
    options.fs (1,1) double {mustBePositive} = 51200
    options.beta (1,1) double {mustBePositive} = 1/4
    options.gamma (1,1) double {mustBePositive} = 1/2
end

if isfield(p,'contact') && ~strcmpi(p.contact.model,'linear')
    error('pg_simulate_loaded_v2_response:NonlinearContact', ...
        ['The linear loaded solver cannot apply contact loss consistently. ' ...
         'Use pg_simulate_nonlinear_contact_response instead.']);
end
if isfield(p,'contact') && p.contact.friction.enabled
    error('pg_simulate_loaded_v2_response:FrictionNotCalibrated', ...
        'Friction is disabled until its parameters are calibrated.');
end

fs = options.fs;
dt = 1/fs;
nSample = floor(duration*fs)+1;
time = (0:nSample-1)'*dt;
nSecondOrder = p.map.n+numel(model.frequencyHz);
Cacc = pg_sensor_observation(p,model);
sensor = pg_sensor_configuration(p);

[M,C,K,meta] = pg_assemble_coupled(time(1),p,fault,model);
[force,forceDetail] = pg_coupled_force_loaded_v2(time(1),p,meta);
q = K\force;
velocity = zeros(nSecondOrder,1);
acceleration = M\(force-C*velocity-K*q);

sensorAcceleration = zeros(nSample,sensor.nSensor);
meshForceSp = zeros(nSample,p.model.nPlanet);
meshForceRp = zeros(nSample,p.model.nPlanet);
inputTorque = zeros(nSample,1);
faultLossSunPlanet = zeros(nSample,p.model.nPlanet);
faultLossRingPlanet = zeros(nSample,p.model.nPlanet);
sensorAcceleration(1,:) = (Cacc*acceleration).';
contact = pg_coupled_contact_force(time(1),q,velocity,p,meta);
meshForceSp(1,:) = contact.sunPlanet;
meshForceRp(1,:) = contact.ringPlanet;
inputTorque(1) = forceDetail.operating.inputTorque;
faultLossSunPlanet(1,:) = meta.tvms.faultStiffnessLossSunPlanet;
faultLossRingPlanet(1,:) = meta.tvms.faultStiffnessLossRingPlanet;

beta = options.beta;
gamma = options.gamma;
for k = 2:nSample
    qPredict = q+dt*velocity+dt^2*(0.5-beta)*acceleration;
    vPredict = velocity+dt*(1-gamma)*acceleration;

    [M,C,K,meta] = pg_assemble_coupled(time(k),p,fault,model);
    [force,forceDetail] = pg_coupled_force_loaded_v2(time(k),p,meta);
    effectiveMass = M+gamma*dt*C+beta*dt^2*K;
    accelerationNew = effectiveMass\ ...
        (force-C*vPredict-K*qPredict);

    q = qPredict+beta*dt^2*accelerationNew;
    velocity = vPredict+gamma*dt*accelerationNew;
    acceleration = accelerationNew;
    sensorAcceleration(k,:) = (Cacc*acceleration).';
    contact = pg_coupled_contact_force(time(k),q,velocity,p,meta);
    meshForceSp(k,:) = contact.sunPlanet;
    meshForceRp(k,:) = contact.ringPlanet;
    inputTorque(k) = forceDetail.operating.inputTorque;
    faultLossSunPlanet(k,:) = meta.tvms.faultStiffnessLossSunPlanet;
    faultLossRingPlanet(k,:) = meta.tvms.faultStiffnessLossRingPlanet;
end

motion = pg_motion_state(time,p);
faultTruth = pg_fault_event_truth(time, ...
    faultLossSunPlanet+faultLossRingPlanet,p,fault);

response.time = time;
response.fs = fs;
response.acceleration = sensorAcceleration;
response.sensorAnglesDeg = sensor.locationAnglesDeg;
response.sensorSensitiveAnglesDeg = sensor.sensitiveAnglesDeg;
response.meshForceSunPlanet = meshForceSp;
response.meshForceRingPlanet = meshForceRp;
response.inputTorque = inputTorque;
response.carrierAngleRad = motion.phi_c(:);
response.carrierSpeedRadS = motion.omega_c(:);
response.meshPhaseRad = motion.meshPhaseRad(:);
response.meshFrequencyHz = motion.f_mesh(:);
response.faultStiffnessLossSunPlanet = faultLossSunPlanet;
response.faultStiffnessLossRingPlanet = faultLossRingPlanet;
response.faultEventTruth = faultTruth;
response.outputTorque = p.operating.outputTorque;
response.fault = fault;
response.modelType = sprintf('loaded 18-DOF LPM plus %d elastic-ring modes', ...
    numel(model.frequencyHz));
response.integration = 'Newmark average acceleration';
response.assumption = p.v2.assumption;
if isfield(p,'v3') && isfield(p.v3,'assumption')
    response.assumption = p.v3.assumption;
end
end
