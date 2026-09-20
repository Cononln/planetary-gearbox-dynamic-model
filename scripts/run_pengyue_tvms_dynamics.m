function result = run_pengyue_tvms_dynamics(options)
%RUN_PENGYUE_TVMS_DYNAMICS Run the 18-DOF Peng Yue TVMS dynamics route.
% The legacy route is untouched.  This runner prepares the Peng Yue lookup
% once, carries it in p.mesh.pengyueSpModel, and uses the existing ODE45/FFT
% response definition (sun-planet DTE of the unified 18-DOF system).

arguments
    options.Case (1,1) string {mustBeMember(options.Case,["healthy","crack45"])} = "healthy"
    options.OutputDir (1,1) string = ""
    options.CrackDepthMm (1,1) double {mustBeNonnegative} = 3
    options.FoundationModel (1,1) string = "none"
    options.FoundationComplianceSun (1,1) double {mustBeNonnegative} = 0
    options.FoundationCompliancePlanet (1,1) double {mustBeNonnegative} = 0
    options.FoundationInnerRadiusSun (1,1) double = NaN
    options.FoundationInnerRadiusPlanet (1,1) double = NaN
    options.FoundationRootHalfAngleSun (1,1) double = NaN
    options.FoundationRootHalfAnglePlanet (1,1) double = NaN
    options.FoundationUfMode (1,1) string = "radial"
    options.RecordCarrierCycles (1,1) double {mustBeInteger,mustBePositive} = 6
    options.DiscardCarrierCycles (1,1) double {mustBeInteger,mustBeNonnegative} = 1
    options.RelTol (1,1) double {mustBePositive} = 1e-5
    options.AbsTol (1,1) double {mustBePositive} = 1e-8
    options.MaxStepMultiplier (1,1) double {mustBePositive} = 12
end

setup_paths;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
if strlength(options.OutputDir)==0
    outDir = fullfile(projectRoot,'results','pengyue_tvms_dynamics');
else
    outDir = char(options.OutputDir);
end
if ~isfolder(outDir); mkdir(outDir); end
p0 = pg_parameters_pengyue_fig314;
assert(p0.map.n==18,'This runner requires the 18-DOF model.');

fHealthy = struct('type','none','tooth',1,'crackDepth',0,'crackAngle',0);
if options.Case == "healthy"
    fault = fHealthy;
else
    fault = struct('type','sun','tooth',1,'crackDepth',options.CrackDepthMm*1e-3, ...
        'crackAngle',deg2rad(45),'crackAngleDeg',45);
end
[p,model] = pg_prepare_pengyue_tvms(p0,fault, ...
    'FoundationModel',options.FoundationModel, ...
    'FoundationComplianceSun',options.FoundationComplianceSun, ...
    'FoundationCompliancePlanet',options.FoundationCompliancePlanet, ...
    'FoundationInnerRadiusSun',options.FoundationInnerRadiusSun, ...
    'FoundationInnerRadiusPlanet',options.FoundationInnerRadiusPlanet, ...
    'FoundationRootHalfAngleSun',options.FoundationRootHalfAngleSun, ...
    'FoundationRootHalfAnglePlanet',options.FoundationRootHalfAnglePlanet, ...
    'FoundationUfMode',options.FoundationUfMode);
fm = p.kin.f_mesh;
carrierPeriod = 1/p.kin.f_c;
discardTime = options.DiscardCarrierCycles*carrierPeriod;
recordDuration = options.RecordCarrierCycles*carrierPeriod;
totalDuration = discardTime+recordDuration;
maxStep = 1/(options.MaxStepMultiplier*fm);

fprintf('Peng Yue 18-DOF TVMS dynamics: %s\n',options.Case);
fprintf('  fm=%.12f Hz, fs_fault=%.12f Hz, DOF=%d\n', ...
    fm,p.kin.f_sun_fault,p.map.n);
fprintf('  total=%.6f s, discard=%.6f s, record=%.6f s, MaxStep=%.6g s\n', ...
    totalDuration,discardTime,recordDuration,maxStep);
fprintf('  crackDepth=%.6g mm, foundationModel=%s, c_s=%.6g, c_p=%.6g m/N\n', ...
    options.CrackDepthMm,options.FoundationModel, ...
    options.FoundationComplianceSun,options.FoundationCompliancePlanet);

r = pg_simulate_pengyue_dte(totalDuration,p,fault,'fs',p.model.fs, ...
    'RelTol',options.RelTol,'AbsTol',options.AbsTol,'MaxStep',maxStep);
keep = r.time>=discardTime;
r.time = r.time(keep)-discardTime;
r.q = r.q(keep,:);
r.qd = r.qd(keep,:);
r.dteSunPlanet = r.dteSunPlanet(keep,:);
r.transmissionErrorSunPlanet = r.transmissionErrorSunPlanet(keep,:);

[kSp,detail] = pg_eval_pengyue_sp_tvms(r.time,p,fault,model);
[kRp,detailRp] = pg_eval_pengyue_pr_tvms(r.time,p,fault,p.mesh.pengyuePrModel);
x = r.dteSunPlanet(:,1);
x = x-mean(x);
n = numel(x);
freq = (0:floor(n/2))'*r.fs/n;
X = fft(x);
amp = 2*abs(X(1:numel(freq)))/n;
amp(1) = abs(X(1))/n;

% Local peak near each theoretical mesh harmonic; no lines are injected.
orders = (1:9).';
theory = orders*fm;
measured = nan(size(theory)); amplitude = nan(size(theory));
for io=1:numel(orders)
    mask = abs(freq-theory(io))<=max(2,3*(freq(2)-freq(1)));
    if any(mask)
        idx = find(mask);
        [amplitude(io),jj] = max(amp(idx));
        measured(io) = freq(idx(jj));
    end
end

stem = char(options.Case);
csvResponse = fullfile(outDir,[stem '_response.csv']);
csvSpectrum = fullfile(outDir,[stem '_spectrum.csv']);
csvOrders = fullfile(outDir,[stem '_mesh_order_metrics.csv']);
T = table(r.time,r.dteSunPlanet(:,1),kSp(:,1),kSp(:,2),kSp(:,3), ...
    kRp(:,1),kRp(:,2),kRp(:,3), ...
    detail.activePairCount(:,1),detail.faultPairCount(:,1), ...
    detailRp.activePairCount(:,1), ...
    'VariableNames',{'time_s','dte_sp1_m','k_sp1_Npm','k_sp2_Npm','k_sp3_Npm', ...
    'k_rp1_Npm','k_rp2_Npm','k_rp3_Npm', ...
    'active_pair_count_sp1','fault_pair_count_sp1','active_pair_count_rp1'});
writetable(T,csvResponse);
writetable(table(freq,amp,'VariableNames',{'frequency_Hz','amplitude_m'}),csvSpectrum);
writetable(table(orders,theory,measured,amplitude, ...
    'VariableNames',{'order','theory_Hz','measured_Hz','amplitude_m'}),csvOrders);

r.kSunPlanet = kSp;
r.kPlanetRing = kRp;
r.tvmsDetail = detail;
r.tvmsDetailPlanetRing = detailRp;
r.frequencyHz = freq;
r.spectrumAmplitude = amp;
r.meshOrderMetrics = table(orders,theory,measured,amplitude, ...
    'VariableNames',{'order','theory_Hz','measured_Hz','amplitude_m'});
r.tvmsSource = ['Peng Yue potential-energy SP and PR TVMS; PR ring-foundation ' ...
    'geometry missing and marked provisional'];
r.tvmsModel = model;
r.steadyStartOriginalTime = discardTime;
r.recordDuration = recordDuration;
r.solverOptions = struct('RelTol',options.RelTol,'AbsTol',options.AbsTol, ...
    'MaxStep',maxStep,'recordCarrierCycles',options.RecordCarrierCycles, ...
    'discardCarrierCycles',options.DiscardCarrierCycles);
save(fullfile(outDir,[stem '_response.mat']),'r','-v7.3');

result = struct('case',stem,'dof',p.map.n,'meshFrequencyHz',fm, ...
    'faultFrequencyHz',p.kin.f_sun_fault,'responseCsv',csvResponse, ...
    'spectrumCsv',csvSpectrum,'orderMetricsCsv',csvOrders, ...
    'meshOrderMetrics',r.meshOrderMetrics);
fprintf('  response: %s\n  spectrum: %s\n  orders: %s\n', ...
    csvResponse,csvSpectrum,csvOrders);
disp(r.meshOrderMetrics);
end
