function bundle = export_pengyue_tvms_for_codex(options)
%EXPORT_PENGYUE_TVMS_FOR_CODEX
% Recompute Peng-Yue sun-planet TVMS and export numerical arrays ONLY.
% No figure is produced. Codex can combine these arrays/functions with the
% existing 18-DOF dynamic solver and regenerate the response spectrum.
%
% Files written:
%   results/pengyue_tvms_bundle.mat
%   results/pengyue_tvms_one_cycle.csv
%
% A single optional GlobalScale is allowed. If Codex calibrates absolute
% stiffness against Peng Yue Fig. 3-2, the SAME scale must be used for
% healthy and all crack angles.

arguments
    options.GlobalScale (1,1) double {mustBePositive} = 1
    options.FoundationModel (1,1) string = "none"
    options.FoundationComplianceSun (1,1) double {mustBeNonnegative} = 0
    options.FoundationCompliancePlanet (1,1) double {mustBeNonnegative} = 0
    options.FoundationInnerRadiusSun (1,1) double = NaN
    options.FoundationInnerRadiusPlanet (1,1) double = NaN
    options.FoundationRootHalfAngleSun (1,1) double = NaN
    options.FoundationRootHalfAnglePlanet (1,1) double = NaN
    options.FoundationUfMode (1,1) string = "radial"
    options.CrackDepthMm (1,1) double {mustBeNonnegative} = 3
    options.SamplesPerMeshCycle (1,1) double {mustBeInteger,mustBeGreaterThan(options.SamplesPerMeshCycle,200)} = 2001
end

setup_paths;
p = pg_parameters_pengyue_fig314;
angles = [15 30 45 60 75];
% The current Fig. 3-2 validation uses q=3 mm.  The later Fig. 3-14 angle
% sweep can still be reproduced explicitly with 'CrackDepthMm',4.
q = options.CrackDepthMm*1e-3;

% One reference mesh cycle. This is sufficient to inspect healthy waveform;
% a long z_s-mesh-period record is also exported to retain fault-tooth events.
u1 = linspace(0,1,options.SamplesPerMeshCycle).';
t1 = u1/p.kin.f_mesh;

% Healthy model.
fHealthy = struct('type','none','tooth',1,'crackDepth',0,'crackAngle',0);
mdlHealthy = pg_build_pengyue_sp_tvms_model(p,fHealthy, ...
    'GlobalScale',options.GlobalScale, ...
    'FoundationModel',options.FoundationModel, ...
    'FoundationComplianceSun',options.FoundationComplianceSun, ...
    'FoundationCompliancePlanet',options.FoundationCompliancePlanet, ...
    'FoundationInnerRadiusSun',options.FoundationInnerRadiusSun, ...
    'FoundationInnerRadiusPlanet',options.FoundationInnerRadiusPlanet, ...
    'FoundationRootHalfAngleSun',options.FoundationRootHalfAngleSun, ...
    'FoundationRootHalfAnglePlanet',options.FoundationRootHalfAnglePlanet, ...
    'FoundationUfMode',options.FoundationUfMode);
[kHealthy,detailHealthy] = pg_eval_pengyue_sp_tvms(t1,p,fHealthy,mdlHealthy);

models = cell(numel(angles),1);
kFaultOneCycle = zeros(numel(t1),numel(angles));
for ia = 1:numel(angles)
    f = struct('type','sun','tooth',1,'crackDepth',q, ...
        'crackAngle',deg2rad(angles(ia)),'crackAngleDeg',angles(ia));
    models{ia} = pg_build_pengyue_sp_tvms_model(p,f, ...
        'GlobalScale',options.GlobalScale, ...
        'FoundationModel',options.FoundationModel, ...
        'FoundationComplianceSun',options.FoundationComplianceSun, ...
        'FoundationCompliancePlanet',options.FoundationCompliancePlanet, ...
        'FoundationInnerRadiusSun',options.FoundationInnerRadiusSun, ...
        'FoundationInnerRadiusPlanet',options.FoundationInnerRadiusPlanet, ...
        'FoundationRootHalfAngleSun',options.FoundationRootHalfAngleSun, ...
        'FoundationRootHalfAnglePlanet',options.FoundationRootHalfAnglePlanet, ...
        'FoundationUfMode',options.FoundationUfMode);

    % Shift the reference time so tooth 1 is active in planet 1. This gives
    % a directly comparable one-cycle waveform for all gamma values.
    [kf,~] = pg_eval_pengyue_sp_tvms(t1,p,f,models{ia});
    kFaultOneCycle(:,ia) = kf(:,1);
end

% Long record covering one sun-tooth repeat for planet 1 (z_s mesh cycles).
nLong = p.gear.zs*500+1;
uLong = linspace(0,p.gear.zs,nLong).';
tLong = uLong/p.kin.f_mesh;
kLong = cell(numel(angles),1);
for ia = 1:numel(angles)
    f = struct('type','sun','tooth',1,'crackDepth',q, ...
        'crackAngle',deg2rad(angles(ia)),'crackAngleDeg',angles(ia));
    [kk,dd] = pg_eval_pengyue_sp_tvms(tLong,p,f,models{ia});
    kLong{ia}.time = tLong;
    kLong{ia}.meshCounter = uLong;
    kLong{ia}.kSunPlanet = kk;
    kLong{ia}.detail = dd;
end

bundle.parameters = p;
bundle.anglesDeg = angles;
bundle.crackDepthM = q;
bundle.healthy.model = mdlHealthy;
bundle.healthy.timeOneCycle = t1;
bundle.healthy.meshCounterOneCycle = u1;
bundle.healthy.kSunPlanetOneCycle = kHealthy;
bundle.healthy.detailOneCycle = detailHealthy;
bundle.faultModels = models;
bundle.kFaultPlanet1OneCycle = kFaultOneCycle;
bundle.longFaultRecords = kLong;
bundle.globalScale = options.GlobalScale;
bundle.note = ['Use kSunPlanet(t) directly as the time-varying stiffness in ' ...
    'the existing 18-DOF assembly. Do not convert it to an arbitrary ' ...
    'sinusoid and do not separately rescale different gamma cases.'];

if ~exist('results','dir'); mkdir('results'); end
save(fullfile('results','pengyue_tvms_bundle.mat'),'bundle','-v7.3');

T = table(u1,t1,kHealthy(:,1), ...
    kFaultOneCycle(:,1),kFaultOneCycle(:,2),kFaultOneCycle(:,3), ...
    kFaultOneCycle(:,4),kFaultOneCycle(:,5), ...
    'VariableNames',{'mesh_cycle','time_s','healthy_Npm', ...
    'gamma15_Npm','gamma30_Npm','gamma45_Npm','gamma60_Npm','gamma75_Npm'});
writetable(T,fullfile('results','pengyue_tvms_one_cycle.csv'));

fprintf('Peng-Yue TVMS numerical export finished. No figure generated.\n');
fprintf('  contact ratio eps_alpha = %.6f\n',mdlHealthy.contactRatio);
fprintf('  raw/global scale        = %.6g\n',options.GlobalScale);
fprintf('  MAT: results/pengyue_tvms_bundle.mat\n');
fprintf('  CSV: results/pengyue_tvms_one_cycle.csv\n');
end
