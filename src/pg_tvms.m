function [stiffness, detail] = pg_tvms(t, p, fault)
%PG_TVMS Time-varying mesh stiffness with tooth-level fault bookkeeping.
% A tooth pair stays active for contactRatio mesh periods. The stiffness of
% each active pair is integrated over retained face-width slices. A broken
% tooth therefore weakens only its own pair, while any overlapping healthy
% pair remains intact.

if nargin < 3
    fault = pg_fault_case('none',0,p);
end

tShape = size(t);
t = t(:);
nTime = numel(t);
nPlanet = p.model.nPlanet;
stiffness.sunPlanet = zeros(nTime,nPlanet);
stiffness.ringPlanet = zeros(nTime,nPlanet);
detail.faultPairCountSunPlanet = zeros(nTime,nPlanet);
detail.faultPairCountRingPlanet = zeros(nTime,nPlanet);
detail.faultStiffnessLossSunPlanet = zeros(nTime,nPlanet);
detail.faultStiffnessLossRingPlanet = zeros(nTime,nPlanet);

motion = pg_motion_state(t,p);
meshCounter = motion.meshCounter(:);
if isfield(fault,'speedRippleFraction') && fault.speedRippleFraction > 0
    % Integral of the reference-paper first-order speed ripple. The
    % subtraction preserves the requested mesh phase at t=0.
    omegaRipple = 2*pi*p.kin.f_s;
    phase = fault.speedRipplePhase;
    meshCounter = meshCounter + fault.speedRippleFraction*p.kin.f_mesh/ ...
        omegaRipple*(sin(omegaRipple*t+phase)-sin(phase));
end

for it = 1:nTime
    u = meshCounter(it);
    for ip = 1:nPlanet
        sunOffset = (ip-1)*(p.gear.zs/nPlanet);
        [ksp,nFaultSp,lossSp] = oneMesh(u, ...
            p.tvms.sunPlanet.contactRatio,p.mesh.sunPlanet.kMean, ...
            p.gear.zs,p.gear.zp,sunOffset,0,ip,'sunPlanet',p,fault);
        ringOffset = (ip-1)*(p.gear.zr/nPlanet);
        [krp,nFaultRp,lossRp] = oneMesh(u, ...
            p.tvms.ringPlanet.contactRatio,p.mesh.ringPlanet.kMean, ...
            p.gear.zr,p.gear.zp,ringOffset, ...
            p.tvms.planetRingToothOffset,ip,'ringPlanet',p,fault);
        stiffness.sunPlanet(it,ip) = ksp;
        stiffness.ringPlanet(it,ip) = krp;
        detail.faultPairCountSunPlanet(it,ip) = nFaultSp;
        detail.faultPairCountRingPlanet(it,ip) = nFaultRp;
        detail.faultStiffnessLossSunPlanet(it,ip) = lossSp;
        detail.faultStiffnessLossRingPlanet(it,ip) = lossRp;
    end
end

if isscalar(tShape) || prod(tShape)==1
    stiffness.sunPlanet = reshape(stiffness.sunPlanet,1,nPlanet);
    stiffness.ringPlanet = reshape(stiffness.ringPlanet,1,nPlanet);
    detail.faultPairCountSunPlanet = ...
        reshape(detail.faultPairCountSunPlanet,1,nPlanet);
    detail.faultPairCountRingPlanet = ...
        reshape(detail.faultPairCountRingPlanet,1,nPlanet);
    detail.faultStiffnessLossSunPlanet = ...
        reshape(detail.faultStiffnessLossSunPlanet,1,nPlanet);
    detail.faultStiffnessLossRingPlanet = ...
        reshape(detail.faultStiffnessLossRingPlanet,1,nPlanet);
end
end

function [kTotal,nFault,stiffnessLoss] = oneMesh(u,epsilon,kMean,zCenter,zPlanet, ...
    centerOffset,planetOffset,planetIndex,meshType,p,fault)

nEntry = floor(u);
xi = u-nEntry;
pairEntry = nEntry;
pairAge = xi;
if xi < epsilon-1
    pairEntry = [pairEntry,nEntry-1];
    pairAge = [pairAge,1+xi];
end

% The mean of all active profile contributions over one mesh period is
% epsilon*(profileFloor + profileSin2/2).
meanProfileSum = epsilon*(p.tvms.profileFloor + p.tvms.profileSin2/2);
kPairReference = kMean/meanProfileSum;
kTotal = 0;
nFault = 0;
stiffnessLoss = 0;

for jp = 1:numel(pairEntry)
    entry = pairEntry(jp);
    age = pairAge(jp);
    shape = p.tvms.profileFloor + p.tvms.profileSin2* ...
        sin(pi*age/epsilon)^2;

    centerTooth = mod(entry+centerOffset,zCenter)+1;
    planetTooth = mod(-entry+planetOffset,zPlanet)+1;
    isFaultPair = false;
    if strcmp(fault.type,'sun') && strcmp(meshType,'sunPlanet')
        isFaultPair = centerTooth == fault.tooth;
    elseif strcmp(fault.type,'planet') && planetIndex == fault.targetPlanet
        isFaultPair = planetTooth == fault.tooth;
    end

    retainedRatio = faultPairStiffnessRatio(isFaultPair,age/epsilon,p,fault);
    kTotal = kTotal + kPairReference*shape*retainedRatio;
    nFault = nFault + double(isFaultPair);
    stiffnessLoss = stiffnessLoss + ...
        kPairReference*shape*(1-retainedRatio);
end
end

function ratio = faultPairStiffnessRatio(isFaultPair,contactFraction,p,fault)
damageModel = 'partial_width_break';
if isfield(fault,'damageModel')
    damageModel = fault.damageModel;
end
if ~isFaultPair || strcmp(damageModel,'none')
    ratio = 1;
    return;
end

if strcmp(damageModel,'root_crack_energy')
    ratio = interp1(fault.energyLookup.contactFraction, ...
        fault.energyLookup.stiffnessRatio,contactFraction,'linear','extrap');
    ratio = min(max(ratio,0),1);
    return;
end

if ~strcmp(damageModel,'partial_width_break')
    error('pg_tvms:DamageModel','Unsupported damage model: %s',damageModel);
end

% Midpoint integration of the local tooth-pair stiffness across face width.
b = p.gear.faceWidthSunPlanet;
nSlice = p.tvms.nWidthSlices;
z = ((1:nSlice)-0.5)*(b/nSlice);
localComplianceInverse = ones(size(z));
healthyIntegral = sum(localComplianceInverse)*(b/nSlice);
if fault.severity > 0
    retained = z >= fault.missingFaceWidth;
else
    retained = true(size(z));
end
remainingIntegral = sum(localComplianceInverse(retained))*(b/nSlice);
ratio = remainingIntegral/healthyIntegral;
end
