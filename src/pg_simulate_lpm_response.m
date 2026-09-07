function response = pg_simulate_lpm_response(duration,p,fault,options)
%PG_SIMULATE_LPM_RESPONSE Integrate the default 18-DOF LPM sensor model.
% Gear dynamics use Newmark average acceleration.  Two one-way local
% mass-spring-damper observers convert path-weighted ring-mesh forces to the
% fixed 0/90-deg accelerometer responses without an FE/modal ring model.

arguments
    duration (1,1) double {mustBePositive}
    p struct
    fault struct
    options.fs (1,1) double {mustBePositive} = 51200
    options.beta (1,1) double {mustBePositive} = 1/4
    options.gamma (1,1) double {mustBePositive} = 1/2
end

if isfield(p,'contact') && ~strcmpi(p.contact.model,'linear')
    error('pg_simulate_lpm_response:NonlinearContact', ...
        ['The linear LPM solver cannot apply contact loss consistently. ' ...
         'Use pg_simulate_nonlinear_contact_response instead.']);
end
if isfield(p,'contact') && p.contact.friction.enabled
    error('pg_simulate_lpm_response:FrictionNotCalibrated', ...
        'Friction is disabled until its parameters are calibrated.');
end

fs = options.fs;
dt = 1/fs;
nSample = floor(duration*fs)+1;
time = (0:nSample-1)'*dt;
nq = p.map.n;
sensor = pg_sensor_configuration(p);
nSensor = sensor.nSensor;
CaccRigid = pg_lpm_sensor_observation(p);

q = zeros(nq,1);
velocity = zeros(nq,1);
pathQ = zeros(nSensor,1);
pathV = zeros(nSensor,1);
sensorAcceleration = zeros(nSample,nSensor);
meshForceSp = zeros(nSample,p.model.nPlanet);
meshForceRp = zeros(nSample,p.model.nPlanet);
pathWeights = zeros(nSample,nSensor,p.model.nPlanet);
inputTorque = zeros(nSample,1);

pathMass = p.sensor.pathEffectiveMass;
pathOmega = 2*pi*p.sensor.pathNaturalFrequencyHz;
pathDamping = 2*p.sensor.pathZeta*pathOmega*pathMass;
pathStiffness = pathMass*pathOmega^2;

[M,C,K,meta] = pg_assemble_system(time(1),p,fault);
[teForce,~] = pg_baseline_force(time(1),meta.phi_c,p,meta);
[operatingForce,operatingDetail] = pg_operating_force(time(1),p,nq);
externalForce = teForce+operatingForce;
isLoaded = abs(operatingDetail.nominalInputTorque)>0 || ...
    abs(operatingDetail.outputTorque)>0;
if isLoaded
    % Establish the static loaded mesh deflection. Subsequent TVMS changes
    % act on this nonzero preload and generate the physical crack force
    % equivalent to -DeltaK(t)*q0.
    q = K\externalForce;
end
acceleration = M\(externalForce-C*velocity-K*q);
contact = pg_mesh_contact_force(time(1),q,velocity,p,meta);
[pathForce,weights] = pg_lpm_path_force(time(1),contact,p);
if isLoaded
    pathQ = pathForce/pathStiffness;
    pathA = zeros(nSensor,1);
else
    pathA = (pathForce-pathDamping*pathV-pathStiffness*pathQ)/pathMass;
end
sensorAcceleration(1,:) = (p.sensor.rigidRingContribution* ...
    (CaccRigid*acceleration)+pathA).';
meshForceSp(1,:) = contact.sunPlanet;
meshForceRp(1,:) = contact.ringPlanet;
pathWeights(1,:,:) = reshape(weights,1,nSensor,p.model.nPlanet);
inputTorque(1) = operatingDetail.inputTorque;

beta = options.beta;
gamma = options.gamma;
pathEffectiveMass = pathMass+gamma*dt*pathDamping+ ...
    beta*dt^2*pathStiffness;

for k = 2:nSample
    qPredict = q+dt*velocity+dt^2*(0.5-beta)*acceleration;
    vPredict = velocity+dt*(1-gamma)*acceleration;

    [M,C,K,meta] = pg_assemble_system(time(k),p,fault);
    [teForce,~] = pg_baseline_force(time(k),meta.phi_c,p,meta);
    [operatingForce,operatingDetail] = pg_operating_force(time(k),p,nq);
    externalForce = teForce+operatingForce;
    effectiveMass = M+gamma*dt*C+beta*dt^2*K;
    accelerationNew = effectiveMass\ ...
        (externalForce-C*vPredict-K*qPredict);
    q = qPredict+beta*dt^2*accelerationNew;
    velocity = vPredict+gamma*dt*accelerationNew;
    acceleration = accelerationNew;

    contact = pg_mesh_contact_force(time(k),q,velocity,p,meta);
    [pathForce,weights] = pg_lpm_path_force(time(k),contact,p);
    pathPredict = pathQ+dt*pathV+dt^2*(0.5-beta)*pathA;
    pathVPredict = pathV+dt*(1-gamma)*pathA;
    pathANew = (pathForce-pathDamping*pathVPredict- ...
        pathStiffness*pathPredict)/pathEffectiveMass;
    pathQ = pathPredict+beta*dt^2*pathANew;
    pathV = pathVPredict+gamma*dt*pathANew;
    pathA = pathANew;

    sensorAcceleration(k,:) = (p.sensor.rigidRingContribution* ...
        (CaccRigid*acceleration)+pathA).';
    meshForceSp(k,:) = contact.sunPlanet;
    meshForceRp(k,:) = contact.ringPlanet;
    pathWeights(k,:,:) = reshape(weights,1,nSensor,p.model.nPlanet);
    inputTorque(k) = operatingDetail.inputTorque;
end

response.time = time;
response.fs = fs;
response.acceleration = sensorAcceleration;
response.sensorAnglesDeg = sensor.locationAnglesDeg;
response.sensorSensitiveAnglesDeg = sensor.sensitiveAnglesDeg;
response.meshForceSunPlanet = meshForceSp;
response.meshForceRingPlanet = meshForceRp;
response.pathWeights = pathWeights;
response.inputTorque = inputTorque;
response.outputTorque = p.operating.outputTorque;
response.fault = fault;
if isLoaded
    response.modelType = ['loaded 18-DOF 2-D lumped-parameter gearbox ' ...
        'with local SDOF path observers'];
else
    response.modelType = ['18-DOF 2-D lumped-parameter gearbox with ' ...
        'local SDOF path observers'];
end
response.integration = 'Newmark average acceleration';
response.beta = beta;
response.gamma = gamma;
end
