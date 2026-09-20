function report = validate_pengyue_full_tvms(options)
%VALIDATE_PENGYUE_FULL_TVMS Validate physical SP/PR TVMS before dynamics.
% This script intentionally runs the healthy TVMS only.  It verifies the
% contact ratios, active-pair fractions, Parker phase offsets, and the six
% per-planet stiffness traces.  No crack scaling, filtering, or spectrum
% adjustment is performed.

arguments
    options.OutputDir (1,1) string = ""
    options.Samples (1,1) double {mustBeInteger,mustBeGreaterThan(options.Samples,1000)} = 6001
end

% Resolve all dependencies from this repository checkout.
thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisDir);
if ~isfile(fullfile(projectRoot,'setup_paths.m'))
    error('validate_pengyue_full_tvms:RepositoryRootNotFound', ...
        'Cannot locate setup_paths.m relative to this script.');
end
addpath(projectRoot,fullfile(projectRoot,'src'),'-begin');
if strlength(options.OutputDir)==0
    outDir = fullfile(projectRoot,'results','pengyue_full_tvms_validation');
else
    outDir = char(options.OutputDir);
end
if ~isfolder(outDir); mkdir(outDir); end

p0 = pg_parameters_pengyue_fig314;
assert(p0.map.n==18,'Expected the 18-DOF model.');
fault = struct('type','none','tooth',1,'targetPlanet',1, ...
    'crackDepth',0,'crackAngle',0);
[p,spModel] = pg_prepare_pengyue_tvms(p0,fault, ...
    'FoundationModel','none','LookupPoints',801,'IntegrationPoints',800);
prModel = p.mesh.pengyuePrModel;

u = linspace(0,3,options.Samples).';
t = u/p.kin.f_mesh;
[kSp,dSp] = pg_eval_pengyue_sp_tvms(t,p,fault,spModel);
[kRp,dRp] = pg_eval_pengyue_pr_tvms(t,p,fault,prModel);

% Theoretical double-pair fraction is epsilon-1 for 1<epsilon<2.
spDouble = mean(dSp.activePairCount==2,1);
prDouble = mean(dRp.activePairCount==2,1);

% Phase regression: planet i must equal planet 1 evaluated at the Parker
% local phase offset, modulo one mesh cycle.
spPhase = mod(spModel.meshPhaseCycles,1);
prPhase = mod(prModel.meshPhaseCycles,1);
spPhaseErr = zeros(1,p.model.nPlanet);
prPhaseErr = zeros(1,p.model.nPlanet);
for ip=1:p.model.nPlanet
    if ip==1; continue; end
    [ksShift,~] = pg_eval_pengyue_sp_tvms(t+spPhase(ip)/p.kin.f_mesh,p,fault,spModel);
    [krShift,~] = pg_eval_pengyue_pr_tvms(t+prPhase(ip)/p.kin.f_mesh,p,fault,prModel);
    spPhaseErr(ip) = max(abs(kSp(:,ip)-ksShift(:,1)))/max(max(kSp(:,ip)),eps);
    prPhaseErr(ip) = max(abs(kRp(:,ip)-krShift(:,1)))/max(max(kRp(:,ip)),eps);
end

T = table(t,u,kSp(:,1),kSp(:,2),kSp(:,3),kRp(:,1),kRp(:,2),kRp(:,3), ...
    dSp.activePairCount(:,1),dSp.activePairCount(:,2),dSp.activePairCount(:,3), ...
    dRp.activePairCount(:,1),dRp.activePairCount(:,2),dRp.activePairCount(:,3), ...
    'VariableNames',{'time_s','mesh_phase_cycles','k_sp1_Npm','k_sp2_Npm', ...
    'k_sp3_Npm','k_rp1_Npm','k_rp2_Npm','k_rp3_Npm','active_sp1', ...
    'active_sp2','active_sp3','active_rp1','active_rp2','active_rp3'});
csvPath = fullfile(outDir,'pengyue_full_tvms_healthy.csv');
writetable(T,csvPath);

stats = [min(kSp,[],1); max(kSp,[],1); mean(kSp,1); ...
    min(kRp,[],1); max(kRp,[],1); mean(kRp,1)];
statsNames = {'SP_min';'SP_max';'SP_mean';'PR_min';'PR_max';'PR_mean'};
writetable(array2table(stats,'VariableNames',{'P1','P2','P3'}, ...
    'RowNames',statsNames),fullfile(outDir,'pengyue_full_tvms_stats.csv'), ...
    'WriteRowNames',true);

fid=fopen(fullfile(outDir,'pengyue_full_tvms_report.txt'),'w');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'Peng Yue physical TVMS healthy validation\n');
fprintf(fid,'DOF=%d; fm=%.12g Hz; foundationModel=none\n',p.map.n,p.kin.f_mesh);
fprintf(fid,'SP epsilon=%.12g; PR epsilon=%.12g\n',spModel.contactRatio,prModel.contactRatio);
fprintf(fid,'SP phase cycles=[%.12g %.12g %.12g]\n',spPhase);
fprintf(fid,'PR phase cycles=[%.12g %.12g %.12g]\n',prPhase);
fprintf(fid,'SP full Parker cycles=[%.12g %.12g %.12g]\n',spModel.toothPhaseCycles);
fprintf(fid,'PR full Parker cycles=[%.12g %.12g %.12g]\n',prModel.toothPhaseCycles);
fprintf(fid,'SP double fractions=[%.9g %.9g %.9g], theory=%.9g\n',spDouble,spModel.contactRatio-1);
fprintf(fid,'PR double fractions=[%.9g %.9g %.9g], theory=%.9g\n',prDouble,prModel.contactRatio-1);
fprintf(fid,'SP phase regression relative errors=[%.3g %.3g %.3g]\n',spPhaseErr);
fprintf(fid,'PR phase regression relative errors=[%.3g %.3g %.3g]\n',prPhaseErr);
fprintf(fid,'Ring foundation status: %s\n',prModel.ringFoundationStatus);
fprintf(fid,'Internal geometry status: %s\n',prModel.internalGeometry.geometryStatus);
fprintf(fid,'MISSING Chaari geometry: ring rim/body, tooth-root fillet; sun/planet r_int and theta_f.\n');
fprintf(fid,'SP min/max/mean P1 = %.9g / %.9g / %.9g N/m\n',min(kSp(:,1)),max(kSp(:,1)),mean(kSp(:,1)));
fprintf(fid,'PR min/max/mean P1 = %.9g / %.9g / %.9g N/m\n',min(kRp(:,1)),max(kRp(:,1)),mean(kRp(:,1)));

report=struct('dof',p.map.n,'meshFrequencyHz',p.kin.f_mesh, ...
    'spContactRatio',spModel.contactRatio,'prContactRatio',prModel.contactRatio, ...
    'spPhaseCycles',spPhase,'prPhaseCycles',prPhase, ...
    'spDoublePairFraction',spDouble,'prDoublePairFraction',prDouble, ...
    'spPhaseRelativeError',spPhaseErr,'prPhaseRelativeError',prPhaseErr, ...
    'csvPath',csvPath,'reportPath',fullfile(outDir,'pengyue_full_tvms_report.txt'));

fprintf('Peng Yue physical TVMS validation (healthy only)\n');
fprintf('  DOF=%d; fm=%.9g Hz\n',p.map.n,p.kin.f_mesh);
fprintf('  epsilon_SP=%.9f; epsilon_PR=%.9f\n',spModel.contactRatio,prModel.contactRatio);
fprintf('  SP phases=[%.6f %.6f %.6f]; PR phases=[%.6f %.6f %.6f]\n',spPhase,prPhase);
fprintf('  SP P1 min/max/mean=%.6g/%.6g/%.6g N/m\n',min(kSp(:,1)),max(kSp(:,1)),mean(kSp(:,1)));
fprintf('  PR P1 min/max/mean=%.6g/%.6g/%.6g N/m\n',min(kRp(:,1)),max(kRp(:,1)),mean(kRp(:,1)));
fprintf('  phase errors SP=[%.3g %.3g %.3g], PR=[%.3g %.3g %.3g]\n',spPhaseErr,prPhaseErr);
fprintf('  output=%s\n',outDir);
end
