function report = run_tacholess_joint_phase_demo
%RUN_TACHOLESS_JOINT_PHASE_DEMO Blind dual-channel phase-separation demo.
% Encoder/kinematic phase is used to generate and score the benchmark, but
% it is never passed to the estimators.

setup_paths;
rng(20260829,'twister');
p = pg_parameters_loaded_v2;
model = pg_ring_modal_model(p);
fs = p.model.fs;
durationSec = 5.0;
t = (0:round(durationSec*fs)-1)'/fs;
harmonicOrder = [1,11,18];
centerFrequencyHz = harmonicOrder*p.kin.f_mesh;
halfBandwidthHz = [48,190,280];
harmonicAmplitude = [1.00,0.62,0.42];

% Slowly varying speed is deliberately not carrier synchronous.  This is
% the identifiability condition used to separate it from three-planet path
% orders; the exact phase below remains hidden from every estimator.
carrierFrequencyHz = p.kin.f_c*(1+0.030*sin(2*pi*0.37*t+0.2)+ ...
    0.015*sin(2*pi*1.10*t-0.5)+0.008*sin(2*pi*0.08*t+1.1));
meshFrequencyHz = p.gear.zr*carrierFrequencyHz;
trueMeshPhase = cumtrapz(t,2*pi*meshFrequencyHz);
trueCarrierPhase = trueMeshPhase/p.gear.zr;

% Compute the three-planet coherent moving paths on an angle grid, then
% sample them along the nonuniform true carrier phase.
pathAngle = linspace(0,2*pi,721);
pathField = pg_path_field(centerFrequencyHz,pathAngle,p,model,ones(1,3));
nSample = numel(t);
nSensor = 2;
nHarmonic = numel(harmonicOrder);
truePath = complex(zeros(nSample,nSensor,nHarmonic));
trueComponent = complex(zeros(size(truePath)));
wrappedCarrier = mod(trueCarrierPhase,2*pi);
for ih = 1:nHarmonic
    gridPath = squeeze(pathField.HCombined(ih,:,:));
    pathScale = sqrt(mean(abs(gridPath).^2,'all'));
    gridPath = gridPath/max(pathScale,eps);
    for is = 1:nSensor
        truePath(:,is,ih) = interp1(pathAngle,gridPath(:,is), ...
            wrappedCarrier,'pchip');
        trueComponent(:,is,ih) = harmonicAmplitude(ih)* ...
            truePath(:,is,ih).*exp(1i*harmonicOrder(ih)*trueMeshPhase);
    end
end

cleanSignal = squeeze(real(sum(trueComponent,3)));
snrDb = 14;
noise = zeros(size(cleanSignal));
for is = 1:nSensor
    signalRms = sqrt(mean(cleanSignal(:,is).^2));
    noise(:,is) = signalRms/10^(snrDb/20)*randn(nSample,1);
end
rawSignal = cleanSignal+noise;

% This is the only object delivered to the blind phase method.
extracted = pg_extract_analytic_harmonics(rawSignal,fs, ...
    centerFrequencyHz,halfBandwidthHz);
edgeSec = 0.35;
use = t>=edgeSec & t<=durationSec-edgeSec;
tUse = t(use)-t(find(use,1));
component = extracted.component(use,:,:);
truth.meshPhase = trueMeshPhase(use);
truth.carrierPhase = trueCarrierPhase(use);
truth.meshFrequencyHz = meshFrequencyHz(use);
truth.path = truePath(use,:,:);
truth.encoderCount = mod(floor(truth.carrierPhase/(2*pi)* ...
    p.operating.encoderPpr),p.operating.encoderPpr);

jointOptions.nIteration = 3;
jointOptions.spatialMultiple = p.model.nPlanet;
jointOptions.pathOrder = 4;
jointOptions.ridge = 2e-3;
jointOptions.relaxation = 0.35;
jointOptions.commonAnchorWeight = 0.80;
jointOptions.phaseSmoothSec = 0.012;
jointOptions.frequencySmoothSec = 0.040;
jointOptions.removeCommonSpatialOrder = false;

% Four baselines/ablations and the proposed constrained estimator.
singlePhase = commonPhaseNoPath(component(:,1,1),fs,1, ...
    p.kin.f_mesh,0.012);
multiPhase = commonPhaseNoPath(component,fs,harmonicOrder, ...
    p.kin.f_mesh,0.012);
[relativePhase,relativeCorrected,relativePath] = relativePathOnly( ...
    component,fs,harmonicOrder,p.kin.f_mesh,p.gear.zr,jointOptions);

unconstrainedOptions = jointOptions;
unconstrainedOptions.spatialMultiple = 1;
unconstrainedOptions.pathOrder = 12;
unconstrained = pg_joint_tacholess_phase(component,fs,harmonicOrder, ...
    p.kin.f_mesh,p.gear.zr,unconstrainedOptions);
joint = pg_joint_tacholess_phase(component,fs,harmonicOrder, ...
    p.kin.f_mesh,p.gear.zr,jointOptions);

methodName = {'single_H1';'multiharmonic_no_path'; ...
    'dual_relative_path_only';'unconstrained_spatial_path'; ...
    'three_planet_joint'};
phaseSet = {singlePhase,multiPhase,relativePhase, ...
    unconstrained.meshPhase,joint.meshPhase};
correctedSet = {component,component,relativeCorrected, ...
    unconstrained.correctedComponent,joint.correctedComponent};
pathSet = {[],[],relativePath,unconstrained.pathPhasor,joint.pathPhasor};
parameterCount = [0;0;nHarmonic*(2*jointOptions.pathOrder+1); ...
    nHarmonic*(2*unconstrainedOptions.pathOrder+1); ...
    nHarmonic*(2*jointOptions.pathOrder+1)];

nMethod = numel(methodName);
phaseRmseDeg = zeros(nMethod,1);
phaseP95Deg = zeros(nMethod,1);
meshFrequencyRmseHz = zeros(nMethod,1);
carrierSpeedRmseRpm = zeros(nMethod,1);
dualResidualDeg = zeros(nMethod,1);
phaseConcentration = zeros(nMethod,1);
dualPhaseCoherence = zeros(nMethod,1);
amplitudePreservationError = zeros(nMethod,1);
pathDifferenceRmseDeg = nan(nMethod,1);
for im = 1:nMethod
    [phaseRmseDeg(im),phaseP95Deg(im),meshFrequencyRmseHz(im), ...
        carrierSpeedRmseRpm(im)] = phaseScore(phaseSet{im},truth,fs, ...
        p.gear.zr,0.040);
    dualResidualDeg(im) = dualResidual(correctedSet{im});
    phaseConcentration(im) = componentConcentration(correctedSet{im}, ...
        phaseSet{im},harmonicOrder);
    dualPhaseCoherence(im) = dualCoherence(correctedSet{im});
    amplitudePreservationError(im) = max(abs(abs(correctedSet{im}(:))- ...
        abs(component(:))))/max(abs(component(:)));
    if ~isempty(pathSet{im})
        pathDifferenceRmseDeg(im) = pathDifferenceScore( ...
            pathSet{im},truth.path,component);
    end
end

metrics = table(methodName,parameterCount,phaseRmseDeg,phaseP95Deg, ...
    meshFrequencyRmseHz,carrierSpeedRmseRpm,pathDifferenceRmseDeg, ...
    dualResidualDeg,dualPhaseCoherence,phaseConcentration, ...
    amplitudePreservationError, ...
    'VariableNames',{'method','complex_path_parameter_count', ...
    'mesh_phase_rmse_deg','mesh_phase_p95_deg', ...
    'mesh_frequency_rmse_hz','carrier_speed_rmse_rpm', ...
    'relative_path_phase_rmse_deg','dual_phase_residual_deg', ...
    'dual_phase_coherence','common_phase_concentration', ...
    'amplitude_preservation_error'});
writetable(metrics,fullfile('results','tacholess_joint_phase_metrics.csv'));

% Generic carrier-cycle matrices are saved as the downstream TSA/TPSVD
% interface.  Their amplitudes are untouched by the phase compensation.
rawCarrierMatrix = pg_phase_synchronous_matrix(rawSignal(use,:), ...
    joint.carrierPhase,2048);
correctedCarrierMatrix = pg_phase_synchronous_matrix( ...
    joint.correctedSignal,joint.carrierPhase,2048);

fig = figure('Visible','off','Color','w','Position',[60 50 1220 830]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
plot(tUse,truth.meshFrequencyHz,'k','LineWidth',1.4); hold on;
plot(tUse,instantaneousFrequency(multiPhase,fs,0.040), ...
    'Color',[0.35 0.55 0.85],'LineWidth',0.9);
plot(tUse,joint.instantaneousMeshFrequencyHz, ...
    'Color',[0.85 0.20 0.10],'LineWidth',1.0);
xlim([tUse(1),tUse(end)]); grid on; box off;
xlabel('Time (s)'); ylabel('Mesh frequency (Hz)');
title('Blind instantaneous-speed recovery');
legend('Hidden truth','Multiharmonic, no path','Three-planet joint', ...
    'Location','best');

nexttile;
color = lines(nMethod);
displayMethod = [2,4,5];
for im = displayMethod
    errorDeg = alignedPhaseErrorDeg(phaseSet{im},truth.meshPhase);
    plot(tUse,errorDeg,'Color',color(im,:),'LineWidth',0.8); hold on;
end
xlim([tUse(1),tUse(end)]); ylim([-60 60]); grid on; box off;
xlabel('Time (s)'); ylabel('Mesh-phase error (deg)');
title('Phase error after removal of one constant');
legend(methodName(displayMethod),'Interpreter','none','Location','best');

nexttile;
ihPlot = 2;
truthDifference = angle(truth.path(:,2,ihPlot).* ...
    conj(truth.path(:,1,ihPlot)));
estimatedDifference = angle(joint.pathPhasor(:,2,ihPlot).* ...
    conj(joint.pathPhasor(:,1,ihPlot)));
plot(tUse,rad2deg(truthDifference),'k','LineWidth',1.2); hold on;
plot(tUse,rad2deg(estimatedDifference),'Color',[0.85 0.20 0.10], ...
    'LineWidth',0.9);
xlim([tUse(1),tUse(end)]); ylim([-180 180]); grid on; box off;
xlabel('Time (s)'); ylabel('H_{90}-H_0 phase (deg)');
title(sprintf('Relative moving-path phase, h=%d',harmonicOrder(ihPlot)));
legend('Hidden path truth','Blind estimate','Location','best');

nexttile;
yyaxis left;
methodCategory = categorical(methodName,methodName,'Ordinal',true);
bar(methodCategory,dualResidualDeg,'FaceColor',[0.25 0.55 0.80]);
ylabel('Dual-channel residual (deg)');
yyaxis right;
plot(methodCategory,dualPhaseCoherence,'ro-', ...
    'LineWidth',1.2,'MarkerFaceColor','r');
ylabel('Dual-channel phase coherence'); ylim([0 1.05]);
grid on; box off; title('Differential-path alignment');
exportgraphics(fig,fullfile('results','tacholess_joint_phase_validation.png'), ...
    'Resolution',190);
close(fig);

fig = figure('Visible','off','Color','w','Position',[80 80 1160 650]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;
plot(tUse,rawSignal(use,1),'Color',[0.1 0.35 0.8]); hold on;
plot(tUse,rawSignal(use,2),'Color',[0.85 0.25 0.1]);
xlim([0 min(0.20,tUse(end))]); grid on; box off;
xlabel('Time (s)'); ylabel('Acceleration (normalized)');
title('Blind input: synchronized 0 and 90 degree vibration');
legend('0 deg','90 deg');
nexttile;
plot(tUse,joint.correctedSignal(:,1),'Color',[0.1 0.35 0.8]); hold on;
plot(tUse,joint.correctedSignal(:,2),'Color',[0.85 0.25 0.1]);
xlim([0 min(0.20,tUse(end))]); grid on; box off;
xlabel('Time (s)'); ylabel('Phase-corrected amplitude');
title('Phase-only path compensation (amplitudes retained)');
legend('0 deg','90 deg');
exportgraphics(fig,fullfile('results','tacholess_joint_phase_waveforms.png'), ...
    'Resolution',190);
close(fig);

report.parameters = p;
report.t = tUse;
report.rawSignal = rawSignal(use,:);
report.extractedComponent = component;
report.truth = truth;
report.joint = joint;
report.unconstrained = unconstrained;
report.metrics = metrics;
report.rawCarrierMatrix = rawCarrierMatrix;
report.correctedCarrierMatrix = correctedCarrierMatrix;
report.encoderPolicy = ['encoderCount is retained only in report.truth ', ...
    'for scoring; it is not an input to any phase estimator.'];
save(fullfile('results','tacholess_joint_phase_demo.mat'),'report','-v7.3');

fprintf('Tacholess joint phase demo complete (input SNR %.1f dB).\n',snrDb);
disp(metrics);
end

function phase = commonPhaseNoPath(component,fs,harmonicOrder,f0,smoothSec)
if ismatrix(component)
    component = reshape(component,size(component,1),size(component,2),1);
end
[nSample,nSensor,nHarmonic] = size(component);
harmonicOrder = harmonicOrder(:).';
nominal = 2*pi*f0*(0:nSample-1)'/fs;
numerator = zeros(nSample-1,1);
denominator = zeros(nSample-1,1);
for is = 1:nSensor
    for ih = 1:nHarmonic
        amplitude = abs(component(:,is,ih));
        normalized = amplitude/max(median(amplitude),eps);
        weight = normalized.^2./(1+normalized.^2);
        pairWeight = min(weight(2:end),weight(1:end-1));
        residualIncrement = angle(component(2:end,is,ih).* ...
            conj(component(1:end-1,is,ih)).* ...
            exp(-1i*harmonicOrder(ih)*diff(nominal)))/harmonicOrder(ih);
        numerator = numerator+pairWeight.*residualIncrement;
        denominator = denominator+pairWeight;
    end
end
increment = numerator./max(denominator,eps);
increment = localMovingAverage(increment,round(0.0015*fs));
raw = nominal+[0;cumsum(increment)];
phase = nominal+localMovingAverage(raw-nominal,round(smoothSec*fs));
end

function [phase,corrected,path] = relativePathOnly(component,fs,harmonicOrder, ...
    f0,zr,options)
initial = commonPhaseNoPath(component,fs,harmonicOrder,f0, ...
    options.phaseSmoothSec);
theta = initial/zr;
orders = options.spatialMultiple*(-options.pathOrder:options.pathOrder);
basis = exp(1i*theta*orders);
corrected = component;
path = complex(ones(size(component)));
for ih = 1:numel(harmonicOrder)
    cross = component(:,2,ih).*conj(component(:,1,ih));
    unitCross = cross./max(abs(cross),eps);
    amplitude = sqrt(abs(component(:,1,ih)).*abs(component(:,2,ih)));
    weight = amplitude/max(median(amplitude),eps);
    weight = weight.^2./(1+weight.^2);
    rootWeight = sqrt(weight);
    weightedBasis = basis.*rootWeight;
    normal = weightedBasis'*weightedBasis;
    coefficient = (normal+options.ridge*trace(normal)/size(normal,1)* ...
        eye(size(normal)))\(weightedBasis'*(unitCross.*rootWeight));
    fitted = basis*coefficient;
    fitted = fitted./max(abs(fitted),eps);
    path(:,2,ih) = fitted;
    corrected(:,2,ih) = component(:,2,ih).*conj(fitted);
end
phase = commonPhaseNoPath(corrected,fs,harmonicOrder,f0, ...
    options.phaseSmoothSec);
end

function [rmseDeg,p95Deg,frequencyRmseHz,carrierRmseRpm] = ...
    phaseScore(estimate,truth,fs,zr,smoothSec)
errorDeg = alignedPhaseErrorDeg(estimate,truth.meshPhase);
rmseDeg = sqrt(mean(errorDeg.^2));
p95Deg = percentile(abs(errorDeg),95);
estimatedFrequency = instantaneousFrequency(estimate,fs,smoothSec);
trueFrequency = localMovingAverage(truth.meshFrequencyHz, ...
    round(smoothSec*fs));
edge = round(0.08*fs);
use = (1+edge):(numel(estimate)-edge);
frequencyRmseHz = sqrt(mean((estimatedFrequency(use)-trueFrequency(use)).^2));
carrierRmseRpm = frequencyRmseHz/zr*60;
end

function errorDeg = alignedPhaseErrorDeg(estimate,truth)
difference = estimate-truth;
difference = difference-median(difference);
errorDeg = rad2deg(difference);
end

function frequency = instantaneousFrequency(phase,fs,smoothSec)
frequency = [diff(phase);phase(end)-phase(end-1)]*fs/(2*pi);
frequency = localMovingAverage(frequency,round(smoothSec*fs));
end

function value = dualResidual(component)
[~,~,nHarmonic] = size(component);
residual = zeros(nHarmonic,1);
for ih = 1:nHarmonic
    difference = angle(component(:,2,ih).*conj(component(:,1,ih)));
    center = angle(mean(exp(1i*difference)));
    difference = angle(exp(1i*(difference-center)));
    residual(ih) = sqrt(mean(difference.^2));
end
value = rad2deg(mean(residual));
end

function value = componentConcentration(component,phase,harmonicOrder)
[~,nSensor,nHarmonic] = size(component);
concentration = zeros(nSensor,nHarmonic);
for is = 1:nSensor
    for ih = 1:nHarmonic
        residual = angle(component(:,is,ih))-harmonicOrder(ih)*phase;
        weight = abs(component(:,is,ih));
        concentration(is,ih) = abs(sum(weight.*exp(1i*residual))/sum(weight));
    end
end
value = mean(concentration,'all');
end

function value = dualCoherence(component)
[~,~,nHarmonic] = size(component);
coherence = zeros(nHarmonic,1);
for ih = 1:nHarmonic
    cross = component(:,2,ih).*conj(component(:,1,ih));
    weight = sqrt(abs(component(:,2,ih)).*abs(component(:,1,ih)));
    coherence(ih) = abs(sum(weight.*cross./max(abs(cross),eps))/sum(weight));
end
value = mean(coherence);
end

function rmseDeg = pathDifferenceScore(estimatePath,truePath,component)
[~,~,nHarmonic] = size(component);
sumSquare = 0;
sumWeight = 0;
for ih = 1:nHarmonic
    estimated = angle(estimatePath(:,2,ih).*conj(estimatePath(:,1,ih)));
    truth = angle(truePath(:,2,ih).*conj(truePath(:,1,ih)));
    error = angle(exp(1i*(estimated-truth)));
    weight = sqrt(abs(component(:,1,ih)).*abs(component(:,2,ih)));
    sumSquare = sumSquare+sum(weight.*error.^2);
    sumWeight = sumWeight+sum(weight);
end
rmseDeg = rad2deg(sqrt(sumSquare/sumWeight));
end

function y = localMovingAverage(x,window)
window = max(1,2*floor(max(window,1)/2)+1);
pad = floor(window/2);
xp = [repmat(x(1,:),pad,1);x;repmat(x(end,:),pad,1)];
y = conv2(xp,ones(window,1)/window,'valid');
end

function value = percentile(x,p)
x = sort(x(:));
index = 1+(numel(x)-1)*p/100;
lower = floor(index);
upper = ceil(index);
if lower==upper
    value = x(lower);
else
    value = x(lower)+(index-lower)*(x(upper)-x(lower));
end
end
