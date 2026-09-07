function response = pg_simulate_sensor_response(duration,p,fault,model,options)
%PG_SIMULATE_SENSOR_RESPONSE Integrate the 24-DOF time-varying model.
% Uses the unconditionally stable average-acceleration Newmark method.

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
    error('pg_simulate_sensor_response:NonlinearContact', ...
        ['The linear Newmark solver cannot apply contact loss consistently. ' ...
         'Use pg_simulate_nonlinear_contact_response instead.']);
end
if isfield(p,'contact') && p.contact.friction.enabled
    error('pg_simulate_sensor_response:FrictionNotCalibrated', ...
        'Friction is disabled until its parameters are calibrated.');
end

fs = options.fs;
dt = 1/fs;
nSample = floor(duration*fs)+1;
time = (0:nSample-1)'*dt;
nSecondOrder = p.map.n+numel(model.frequencyHz);
Cacc = pg_sensor_observation(p,model);
sensor = pg_sensor_configuration(p);

q = zeros(nSecondOrder,1);
velocity = zeros(nSecondOrder,1);
acceleration = zeros(nSecondOrder,1);
sensorAcceleration = zeros(nSample,sensor.nSensor);

[M,C,K,meta] = pg_assemble_coupled(time(1),p,fault,model);
force = pg_coupled_force(time(1),p,meta);
acceleration = M\(force-C*velocity-K*q);
sensorAcceleration(1,:) = (Cacc*acceleration).';

beta = options.beta;
gamma = options.gamma;
for k = 2:nSample
    qPredict = q+dt*velocity+dt^2*(0.5-beta)*acceleration;
    vPredict = velocity+dt*(1-gamma)*acceleration;

    [M,C,K,meta] = pg_assemble_coupled(time(k),p,fault,model);
    force = pg_coupled_force(time(k),p,meta);
    effectiveMass = M+gamma*dt*C+beta*dt^2*K;
    accelerationNew = effectiveMass\(force-C*vPredict-K*qPredict);

    q = qPredict+beta*dt^2*accelerationNew;
    velocity = vPredict+gamma*dt*accelerationNew;
    acceleration = accelerationNew;
    sensorAcceleration(k,:) = (Cacc*acceleration).';
end

response.time = time;
response.fs = fs;
response.acceleration = sensorAcceleration;
response.sensorAnglesDeg = sensor.locationAnglesDeg;
response.sensorSensitiveAnglesDeg = sensor.sensitiveAnglesDeg;
response.fault = fault;
response.integration = 'Newmark average acceleration';
response.beta = beta;
response.gamma = gamma;
end
