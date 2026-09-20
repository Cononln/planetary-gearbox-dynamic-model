function [kSp,detail] = pg_eval_pengyue_sp_tvms(t,p,fault,model)
%PG_EVAL_PENGYUE_SP_TVMS Fast evaluation of Peng-Yue sun-planet TVMS.
% Sums one or two simultaneously active tooth-pair stiffnesses according to
% the transverse contact ratio. This is the replacement for the old smooth
% sin^2 TVMS used in pg_tvms.
%
% INPUT
%   t     scalar/vector time [s]
%   p     Peng-Yue parameter struct
%   fault fault struct; type='none' or 'sun'
%   model output of pg_build_pengyue_sp_tvms_model
%
% OUTPUT
%   kSp   [numel(t) x nPlanet] absolute stiffness [N/m]
%
% The local mesh counter for each planet uses Parker's phase relation. A sun
% tooth fault is identified by the integer tooth-pair entry index; therefore
% each individual planet sees the fault once per z_s mesh periods, while the
% three planets are phase shifted and jointly produce N*f_m/z_s.

if nargin < 3 || isempty(fault)
    fault = struct('type','none','tooth',1);
end
if ~isfield(fault,'tooth'); fault.tooth = 1; end

t = t(:);
nT = numel(t);
nP = p.model.nPlanet;
kSp = zeros(nT,nP);
detail.activePairCount = zeros(nT,nP);
detail.faultPairCount = zeros(nT,nP);
detail.localMeshCounter = zeros(nT,nP);
detail.healthyStiffness = zeros(nT,nP);
detail.faultStiffnessLoss = zeros(nT,nP);

epsilon = model.contactRatio;
etaGrid = model.eta;

% Reference mesh counter. If pg_motion_state is available, use it so this
% wrapper stays consistent with the existing repository; otherwise f_m*t.
try
    motion = pg_motion_state(t,p);
    uRef = motion.meshCounter(:);
catch
    uRef = p.kin.f_mesh*t;
end

for ip = 1:nP
    % Unwrapped Parker phase in mesh-cycle units.
    uLocal = uRef + model.toothPhaseCycles(ip);
    detail.localMeshCounter(:,ip) = uLocal;

    for it = 1:nT
        u = uLocal(it);
        nEntry = floor(u);
        xi = u-nEntry;

        entries = nEntry;
        ages = xi;
        if xi < epsilon-1
            entries = [nEntry, nEntry-1]; %#ok<AGROW>
            ages = [xi, 1+xi]; %#ok<AGROW>
        end

        kNow = 0;
        kHealthyNow = 0;
        nf = 0;
        for jp = 1:numel(entries)
            eta = ages(jp)/epsilon;
            eta = min(max(eta,0),1);

            sunTooth = mod(entries(jp),p.gear.zs)+1;
            isFault = strcmpi(fault.type,'sun') && sunTooth==fault.tooth;
            kPairHealthy = interp1(etaGrid,model.kPairHealthy,eta,'linear');
            if isFault
                kPair = interp1(etaGrid,model.kPairFault,eta,'linear');
                nf = nf+1;
            else
                kPair = kPairHealthy;
            end
            kHealthyNow = kHealthyNow+kPairHealthy;
            kNow = kNow+kPair; % Eq. (2-24): active pairs in parallel
        end

        kSp(it,ip) = kNow;
        detail.healthyStiffness(it,ip) = kHealthyNow;
        detail.faultStiffnessLoss(it,ip) = max(kHealthyNow-kNow,0);
        detail.activePairCount(it,ip) = numel(entries);
        detail.faultPairCount(it,ip) = nf;
    end
end
end
