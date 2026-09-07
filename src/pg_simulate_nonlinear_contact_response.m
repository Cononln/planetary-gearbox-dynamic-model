function response = pg_simulate_nonlinear_contact_response( ...
    duration,p,fault,model,options)
%PG_SIMULATE_NONLINEAR_CONTACT_RESPONSE Implicit Newmark contact solution.
% Use this solver for unilateral or backlash contact. Newton iterations use
% the consistent active-contact tangent at every sample.

arguments
    duration (1,1) double {mustBePositive}
    p struct
    fault struct
    model struct
    options.fs (1,1) double {mustBePositive} = 51200
    options.beta (1,1) double {mustBePositive} = 1/4
    options.gamma (1,1) double {mustBePositive} = 1/2
    options.maxNewtonIteration (1,1) double {mustBeInteger,mustBePositive} = 20
    options.relativeTolerance (1,1) double {mustBePositive} = 1e-8
end

if strcmpi(p.contact.model,'linear')
    error('pg_simulate_nonlinear_contact_response:LinearModel', ...
        'Use pg_simulate_loaded_v2_response for the linear contact model.');
end

fs = options.fs;
dt = 1/fs;
time = (0:floor(duration*fs))'*dt;
nSample = numel(time);
n = p.map.n+numel(model.frequencyHz);
sensor = pg_sensor_configuration(p);
Cacc = pg_sensor_observation(p,model);

% A linear loaded equilibrium is a robust initial guess for the nonlinear
% static balance at the first carrier/mesh phase.
[Mlin,Clin,Klin,metaLin] = pg_assemble_coupled(time(1),p,fault,model);
[flin,~] = pg_coupled_force_loaded_v2(time(1),p,metaLin);
q = Klin\flin;
velocity = zeros(n,1);
for iteration = 1:options.maxNewtonIteration
    state = pg_coupled_contact_state(time(1),q,velocity,p,fault,model);
    residual = state.K0*q-state.operatingForce-state.contactForce;
    scale = max(1,norm(state.operatingForce)+norm(state.contactForce));
    if norm(residual) <= options.relativeTolerance*scale
        break;
    end
    q = q-(state.K0+state.contactTangentK)\residual;
end
if iteration == options.maxNewtonIteration && ...
        norm(residual) > options.relativeTolerance*scale
    error('pg_simulate_nonlinear_contact_response:InitialEquilibrium', ...
        'Nonlinear loaded equilibrium did not converge.');
end
state = pg_coupled_contact_state(time(1),q,velocity,p,fault,model);
acceleration = state.M\(state.operatingForce+state.contactForce- ...
    state.C0*velocity-state.K0*q);

sensorAcceleration = zeros(nSample,sensor.nSensor);
meshForceSp = zeros(nSample,p.model.nPlanet);
meshForceRp = zeros(nSample,p.model.nPlanet);
closedSp = false(nSample,p.model.nPlanet);
closedRp = false(nSample,p.model.nPlanet);
newtonIteration = zeros(nSample,1);
sensorAcceleration(1,:) = (Cacc*acceleration).';
meshForceSp(1,:) = state.meshForceSunPlanet;
meshForceRp(1,:) = state.meshForceRingPlanet;
closedSp(1,:) = state.closedSunPlanet;
closedRp(1,:) = state.closedRingPlanet;

beta = options.beta;
gamma = options.gamma;
for it = 2:nSample
    qPredict = q+dt*velocity+dt^2*(0.5-beta)*acceleration;
    vPredict = velocity+dt*(1-gamma)*acceleration;
    accelerationTrial = acceleration;
    for iteration = 1:options.maxNewtonIteration
        qTrial = qPredict+beta*dt^2*accelerationTrial;
        vTrial = vPredict+gamma*dt*accelerationTrial;
        state = pg_coupled_contact_state( ...
            time(it),qTrial,vTrial,p,fault,model);
        residual = state.M*accelerationTrial+state.C0*vTrial+ ...
            state.K0*qTrial-state.operatingForce-state.contactForce;
        scale = max(1,norm(state.operatingForce)+norm(state.contactForce));
        if norm(residual) <= options.relativeTolerance*scale
            break;
        end
        tangent = state.M+gamma*dt*(state.C0+ ...
            state.contactTangentC)+beta*dt^2*(state.K0+ ...
            state.contactTangentK);
        accelerationTrial = accelerationTrial-tangent\residual;
    end
    if iteration == options.maxNewtonIteration && ...
            norm(residual) > options.relativeTolerance*scale
        error('pg_simulate_nonlinear_contact_response:Newton', ...
            'Newton iteration failed at sample %d (t = %.6g s).', ...
            it,time(it));
    end
    acceleration = accelerationTrial;
    q = qPredict+beta*dt^2*acceleration;
    velocity = vPredict+gamma*dt*acceleration;
    newtonIteration(it) = iteration;
    sensorAcceleration(it,:) = (Cacc*acceleration).';
    meshForceSp(it,:) = state.meshForceSunPlanet;
    meshForceRp(it,:) = state.meshForceRingPlanet;
    closedSp(it,:) = state.closedSunPlanet;
    closedRp(it,:) = state.closedRingPlanet;
end

response.time = time;
response.fs = fs;
response.acceleration = sensorAcceleration;
response.sensorAnglesDeg = sensor.locationAnglesDeg;
response.sensorSensitiveAnglesDeg = sensor.sensitiveAnglesDeg;
response.meshForceSunPlanet = meshForceSp;
response.meshForceRingPlanet = meshForceRp;
response.closedSunPlanet = closedSp;
response.closedRingPlanet = closedRp;
response.newtonIteration = newtonIteration;
response.contactModel = p.contact;
response.fault = fault;
response.integration = 'implicit Newmark with active-contact Newton tangent';
end
