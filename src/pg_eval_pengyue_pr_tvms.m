function [kRp,detail] = pg_eval_pengyue_pr_tvms(t,p,fault,model)
%PG_EVAL_PENGYUE_PR_TVMS Evaluate physical internal-mesh TVMS.
% The local PR mesh counter follows Parker gamma_rp^i=Zr*phi_p^i.
% Each active pair is evaluated from the precomputed eta lookup and summed
% directly in parallel. No profileFloor/profileSin2 or extra smoothing is
% used here.

if nargin < 3 || isempty(fault)
    fault = struct('type','none','tooth',1,'crackDepth',0);
end
if ~isfield(fault,'tooth'); fault.tooth = 1; end
t = t(:);
nT = numel(t);
nP = p.model.nPlanet;
kRp = zeros(nT,nP);
detail.activePairCount = zeros(nT,nP);
detail.faultPairCount = zeros(nT,nP);
detail.localMeshCounter = zeros(nT,nP);
detail.healthyStiffness = zeros(nT,nP);
detail.faultStiffnessLoss = zeros(nT,nP);

epsilon = model.contactRatio;
etaGrid = model.eta;
try
    motion = pg_motion_state(t,p);
    uRef = motion.meshCounter(:);
catch
    uRef = p.kin.f_mesh*t;
end

for ip = 1:nP
    % Full Parker phase is retained for tooth identity; fractional phase is
    % what determines the local contact age.  For 17/37/91 this gives
    % [0,1/3,2/3] cycles for the three planets.
    uLocal = uRef + model.toothPhaseCycles(ip);
    detail.localMeshCounter(:,ip) = uLocal;
    for it = 1:nT
        u = uLocal(it);
        nEntry = floor(u);
        xi = u-nEntry;
        entries = nEntry;
        ages = xi;
        if xi < epsilon-1
            entries = [nEntry,nEntry-1]; %#ok<AGROW>
            ages = [xi,1+xi]; %#ok<AGROW>
        end
        kNow = 0;
        kHealthyNow = 0;
        nf = 0;
        for jp = 1:numel(entries)
            eta = min(max(ages(jp)/epsilon,0),1);
            kPairHealthy = interp1(etaGrid,model.kPairHealthy,eta,'linear');
            kPair = kPairHealthy;
            % Ring/planet fault tooth identity is tied to the planet entry
            % index; unlike the legacy fractional ringOffset, this remains
            % integer for every planet and every cycle.
            planetTooth = mod(-entries(jp),p.gear.zp)+1;
            isFault = strcmpi(fault.type,'planet') && ...
                isfield(fault,'targetPlanet') && ip==fault.targetPlanet && ...
                planetTooth==fault.tooth && isfield(fault,'crackDepth') && ...
                fault.crackDepth>0;
            if isFault
                kPair = interp1(etaGrid,model.kPairFault,eta,'linear');
                nf = nf+1;
            end
            kHealthyNow = kHealthyNow+kPairHealthy;
            kNow = kNow+kPair;
        end
        kRp(it,ip) = kNow;
        detail.healthyStiffness(it,ip) = kHealthyNow;
        detail.faultStiffnessLoss(it,ip) = max(kHealthyNow-kNow,0);
        detail.activePairCount(it,ip) = numel(entries);
        detail.faultPairCount(it,ip) = nf;
    end
end
end
