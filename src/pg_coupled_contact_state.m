function state = pg_coupled_contact_state(t,q,qd,p,fault,model)
%PG_COUPLED_CONTACT_STATE Nonlinear mesh forces and consistent tangents.
% Structural supports and ring modes remain linear. TVMS, transmission
% error, and the selected contact law are evaluated at the current state.

motion = pg_motion_state(t,p);
phi_c = motion.phi_c;
if p.contact.friction.enabled
    error('pg_coupled_contact_state:FrictionNotCalibrated', ...
        ['Friction is intentionally disabled until the coefficient and ' ...
         'lubrication condition are calibrated.']);
end
[meshStiffness,tvmsDetail] = pg_tvms(t,p,fault);
zeroStiffness.sunPlanet = zeros(1,p.model.nPlanet);
zeroStiffness.ringPlanet = zeros(1,p.model.nPlanet);
[M,C0,K0,meta] = pg_coupled_matrices( ...
    phi_c,p,model,zeroStiffness);
[operatingForce,operatingDetail] = pg_operating_force(t,p,numel(q));

n = numel(q);
nr = numel(model.frequencyHz);
generalizedContact = zeros(n,1);
Kt = zeros(n);
Ct = zeros(n);
meshForceSp = zeros(1,p.model.nPlanet);
meshForceRp = zeros(1,p.model.nPlanet);
closedSp = false(1,p.model.nPlanet);
closedRp = false(1,p.model.nPlanet);
omegaMesh = motion.omega_mesh;
meshPhaseRad = motion.meshPhaseRad;
dampingOnly = p.contact.dampingOnlyInContact;

for ip = 1:p.model.nPlanet
    gSp = [meta.rigid.mesh(ip).bSunPlanet;zeros(nr,1)];
    gRp = meta.coupling(ip).g;
    kSp = meshStiffness.sunPlanet(ip);
    kRp = meshStiffness.ringPlanet(ip);
    muSp = 1/(gSp.'*(M\gSp));
    muRp = 1/(gRp.'*(M\gRp));
    cSp = 2*p.mesh.sunPlanet.zeta*sqrt(max(kSp,0)*muSp);
    cRp = 2*p.mesh.ringPlanet.zeta*sqrt(max(kRp,0)*muRp);

    phaseSp = p.mesh.sunPlanet.tePhase(ip);
    phaseRp = p.mesh.ringPlanet.tePhase(ip);
    eSp = p.mesh.sunPlanet.teAmplitude*cos(meshPhaseRad+phaseSp);
    eRp = p.mesh.ringPlanet.teAmplitude*cos(meshPhaseRad+phaseRp);
    edSp = -omegaMesh*p.mesh.sunPlanet.teAmplitude* ...
        sin(meshPhaseRad+phaseSp);
    edRp = -omegaMesh*p.mesh.ringPlanet.teAmplitude* ...
        sin(meshPhaseRad+phaseRp);
    % Positive compression is g'*q-e. The legacy linear formulation places
    % g*(k*e+c*edot) on the right and k*g*g' on the left, which is exactly
    % equivalent to the generalized contact force -g*Fcompression.
    deltaSp = gSp.'*q-eSp;
    deltaRp = gRp.'*q-eRp;
    deltaDotSp = gSp.'*qd-edSp;
    deltaDotRp = gRp.'*qd-edRp;

    [fSp,ktSp,ctSp,closedSp(ip)] = pg_contact_law( ...
        deltaSp,deltaDotSp,kSp,cSp,p.contact.model, ...
        p.contact.clearanceSunPlanet,dampingOnly);
    [fRp,ktRp,ctRp,closedRp(ip)] = pg_contact_law( ...
        deltaRp,deltaDotRp,kRp,cRp,p.contact.model, ...
        p.contact.clearanceRingPlanet,dampingOnly);
    generalizedContact = generalizedContact-gSp*fSp-gRp*fRp;
    Kt = Kt+ktSp*(gSp*gSp.')+ktRp*(gRp*gRp.');
    Ct = Ct+ctSp*(gSp*gSp.')+ctRp*(gRp*gRp.');
    meshForceSp(ip) = fSp;
    meshForceRp(ip) = fRp;
end

state.M = M;
state.C0 = C0;
state.K0 = K0;
state.contactForce = generalizedContact;
state.contactTangentK = Kt;
state.contactTangentC = Ct;
state.operatingForce = operatingForce;
state.operatingDetail = operatingDetail;
state.meshForceSunPlanet = meshForceSp;
state.meshForceRingPlanet = meshForceRp;
state.closedSunPlanet = closedSp;
state.closedRingPlanet = closedRp;
state.meta = meta;
state.meta.tvms = tvmsDetail;
end
