function path = pg_apply_fe_ring_transfer_causal(time,meshForceRingPlanet, ...
    carrierAngleRad,p,model,options)
%PG_APPLY_FE_RING_TRANSFER_CAUSAL Causal FE modal state-space path.
%
% The three moving PR forces are projected onto the exported ANSYS modes at
% their instantaneous contact positions.  Assuming SI unit modal mass, every mode is
% then advanced in time as
%
%   eta_ddot + 2*zeta*wn*eta_dot + wn^2*eta = Q_m(t)
%
% using an exact zero-order-hold state transition.  Sensor accelerations
% are recovered from the exported radial mode shapes.  No empirical
% angular-distance weighting, fitted sensor gain, or artificial delay is
% applied.  Differences between sensors therefore arise only from the
% exported modal shapes and the moving contact kinematics.

arguments
    time (:,1) double
    meshForceRingPlanet (:,:) double
    carrierAngleRad (:,1) double
    p struct
    model struct
    options.InitialModalDisplacement double = []
    options.InitialModalVelocity double = []
    options.ReturnModalHistories (1,1) logical = false
end
time = time(:);
carrierAngleRad = carrierAngleRad(:);
n = numel(time);
nMode = numel(model.frequencyHz);
nSensor = size(model.sensorShapeRadial,1);
if n < 2
    error('pg_apply_fe_ring_transfer_causal:Time','At least two samples are required.');
end
if size(meshForceRingPlanet,1)~=n || ...
        size(meshForceRingPlanet,2)~=p.model.nPlanet
    error('pg_apply_fe_ring_transfer_causal:ForceSize','PR force size is invalid.');
end
dt = median(diff(time));
if numel(carrierAngleRad)~=n || any(abs(diff(time)-dt)>1e-8*dt)
    error('pg_apply_fe_ring_transfer_causal:Time','Time/angle vectors are invalid.');
end
if size(model.inputShapeNormal,2)~=nMode || ...
        size(model.sensorShapeRadial,2)~=nMode
    error('pg_apply_fe_ring_transfer_causal:ModeSize','Modal dimensions disagree.');
end

% Moving-contact projection: Q(:,m) is the generalized force of mode m.
Q = zeros(n,nMode);
modalEnergy = zeros(n,nSensor);
sensorShapeEnergy = abs(model.sensorShapeRadial).^2;
contactAngleRad = zeros(n,p.model.nPlanet);
for ip=1:p.model.nPlanet
    theta = pg_fe_ring_contact_angle(carrierAngleRad,p,model,ip);
    contactAngleRad(:,ip)=theta;
    psi=pg_fe_ring_input_shape(theta,model);
    Q = Q + meshForceRingPlanet(:,ip).*psi;
    modalEnergy = modalEnergy + abs(psi).^2*sensorShapeEnergy.';
end

eta = zeros(1,nMode);
etaDot = zeros(1,nMode);
if isempty(options.InitialModalDisplacement)
    eta(:) = 0;
else
    q0=options.InitialModalDisplacement(:).';
    if numel(q0)~=nMode; error('InitialModalDisplacement has wrong size.'); end
    eta=q0;
end
if isempty(options.InitialModalVelocity)
    etaDot(:) = 0;
else
    v0=options.InitialModalVelocity(:).';
    if numel(v0)~=nMode; error('InitialModalVelocity has wrong size.'); end
    etaDot=v0;
end

Ad = zeros(2,2,nMode);
Bd = zeros(2,nMode);
for im=1:nMode
    wn=model.omega(im); zeta=model.zeta(im);
    A=[0 1;-wn^2 -2*zeta*wn]; B=[0;1];
    aug=expm([A B;zeros(1,3)]*dt);
    Ad(:,:,im)=aug(1:2,1:2);
    Bd(:,im)=aug(1:2,3);
end

% Acceleration is evaluated consistently from each modal equation.  Stream
% the modal states so a long high-rate run does not retain three n-by-mode
% history matrices unless explicitly requested.
wn=model.omega(:).'; zeta=model.zeta(:).';
etaDDot = Q(1,:) - 2*zeta.*wn.*etaDot - wn.^2.*eta;
sensorAcceleration=zeros(n,nSensor);
sensorAcceleration(1,:)=etaDDot*model.sensorShapeRadial.';
if options.ReturnModalHistories
    etaHistory=zeros(n,nMode); etaDotHistory=zeros(n,nMode);
    etaDDotHistory=zeros(n,nMode);
    etaHistory(1,:)=eta; etaDotHistory(1,:)=etaDot;
    etaDDotHistory(1,:)=etaDDot;
else
    etaHistory=[]; etaDotHistory=[]; etaDDotHistory=[];
end
ad11=squeeze(Ad(1,1,:)).'; ad12=squeeze(Ad(1,2,:)).';
ad21=squeeze(Ad(2,1,:)).'; ad22=squeeze(Ad(2,2,:)).';
bd1=Bd(1,:); bd2=Bd(2,:);
for k=1:n-1
    etaNew=ad11.*eta+ad12.*etaDot+bd1.*Q(k,:);
    etaDotNew=ad21.*eta+ad22.*etaDot+bd2.*Q(k,:);
    eta=etaNew; etaDot=etaDotNew;
    etaDDot=Q(k+1,:) - 2*zeta.*wn.*etaDot - wn.^2.*eta;
    sensorAcceleration(k+1,:)=etaDDot*model.sensorShapeRadial.';
    if options.ReturnModalHistories
        etaHistory(k+1,:)=eta; etaDotHistory(k+1,:)=etaDot;
        etaDDotHistory(k+1,:)=etaDDot;
    end
end

path.time=time;
path.fs=1/dt;
path.acceleration=sensorAcceleration;
if options.ReturnModalHistories
    path.modalGeneralizedForce=Q;
    path.modalDisplacement=etaHistory;
    path.modalVelocity=etaDotHistory;
    path.modalAcceleration=etaDDotHistory;
end
path.modalEnergy=modalEnergy;
path.contactAngleRad=contactAngleRad;
path.sensorAnglesLocalDeg=model.sensorAnglesLocalDeg;
path.sensorAnglesGlobalDeg=model.sensorAnglesGlobalDeg;
path.sensorObservationStatus=model.sensorObservationStatus;
path.sensorObservationMethod='exported_mode_shapes_causal_state_space';
path.method=['moving PR force projected by exported ANSYS contact mode ' ...
    'shapes; exact-ZOH causal modal state integration; exported radial ' ...
    'sensor mode shapes; no distance weight or fitted gain'];
path.modalHistoriesReturned=options.ReturnModalHistories;
end
