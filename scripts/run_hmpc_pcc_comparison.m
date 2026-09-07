function comparison = run_hmpc_pcc_comparison(healthyFile,faultFile)
%RUN_HMPC_PCC_COMPARISON Compare waveform PCC with phase-only HMPC.
% A 50 s record is band-limited and resampled to 5120 Hz.  For a fault on
% one specified planet, the exact mechanical repeat is 31 carrier turns
% because Zr/Zp = 84/31.  Three processing stages are evaluated:
% nominal-time segmentation, blind common-phase resampling, and common
% phase plus cross-fitted moving-path correction.

setup_paths;
rootDir = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || strlength(string(healthyFile))==0
    healthyFile = fullfile(rootDir,'data_local','bl','healthy_record.MAT');
end
if nargin < 2 || strlength(string(faultFile))==0
    faultFile = fullfile(rootDir,'data_local','planet','pf50_record.MAT');
end

file = {healthyFile;faultFile};
condition = ["BL600_healthy";"PF50_planet_fault"];
stage = ["raw_nominal_phase";"blind_common_phase"; ...
    "common_plus_crossfit_path"];
stageDisplay = ["Nominal time";"Common phase";"Common + path"];
nCondition = numel(file);
nStage = numel(stage);
targetFs = 5120;
segmentStartSec = 5;
segmentDurationSec = 50;
repeatCarrierTurns = 31;
samplesPerRepeat = 32768;
zr = 84;

result = cell(nCondition,1);
for ic = 1:nCondition
    result{ic} = processRecord(file{ic},targetFs,segmentStartSec, ...
        segmentDurationSec,repeatCarrierTurns,samplesPerRepeat,zr);
end

conditionColumn = strings(nCondition*nStage,1);
stageColumn = strings(nCondition*nStage,1);
sensorSubset = strings(nCondition*nStage,1);
blindMeshFrequencyHz = zeros(nCondition*nStage,1);
mechanicalRepeatSec = zeros(nCondition*nStage,1);
repeatCycleCount = zeros(nCondition*nStage,1);
hmpc = zeros(nCondition*nStage,1);
pcc = zeros(nCondition*nStage,1);
cycleCorrelation = zeros(nCondition*nStage,1);
meanCycleStd = zeros(nCondition*nStage,1);
row = 0;
for ic = 1:nCondition
    subsetLabel = strjoin(string(result{ic}.selectedSensorIndex),'+');
    for is = 1:nStage
        row = row+1;
        conditionColumn(row) = condition(ic);
        stageColumn(row) = stage(is);
        sensorSubset(row) = subsetLabel;
        blindMeshFrequencyHz(row) = result{ic}.blindMeshFrequencyHz;
        mechanicalRepeatSec(row) = result{ic}.mechanicalRepeatSec;
        repeatCycleCount(row) = result{ic}.stage(is).cycleCount;
        hmpc(row) = result{ic}.stage(is).hmpc;
        pcc(row) = result{ic}.stage(is).pcc;
        cycleCorrelation(row) = result{ic}.stage(is).cycleCorrelation;
        meanCycleStd(row) = result{ic}.stage(is).meanCycleStd;
    end
end
metrics = table(conditionColumn,stageColumn,sensorSubset, ...
    blindMeshFrequencyHz,mechanicalRepeatSec,repeatCycleCount,hmpc,pcc, ...
    cycleCorrelation,meanCycleStd, ...
    'VariableNames',{'condition','stage','sensor_subset', ...
    'blind_mesh_frequency_hz','mechanical_repeat_sec', ...
    'repeat_cycle_count','hmpc','pcc','cycle_correlation', ...
    'mean_cycle_std'});

hmpcMatrix = reshape(hmpc,nStage,nCondition).';
pccMatrix = reshape(pcc,nStage,nCondition).';
speedGain = normalizedGain(hmpcMatrix(:,1),hmpcMatrix(:,2));
pathGain = normalizedGain(hmpcMatrix(:,2),hmpcMatrix(:,3));
totalGain = normalizedGain(hmpcMatrix(:,1),hmpcMatrix(:,3));
pccSpeedGain = normalizedGain(pccMatrix(:,1),pccMatrix(:,2));
pccPathGain = normalizedGain(pccMatrix(:,2),pccMatrix(:,3));
pccTotalGain = normalizedGain(pccMatrix(:,1),pccMatrix(:,3));
gain = table(condition,speedGain,pathGain,totalGain, ...
    pccSpeedGain,pccPathGain,pccTotalGain, ...
    'VariableNames',{'condition','hmpc_speed_gain', ...
    'hmpc_path_gain','hmpc_total_gain','pcc_speed_gain', ...
    'pcc_path_gain','pcc_total_gain'});
if ~isfolder('results')
    mkdir('results');
end
writetable(metrics,fullfile('results','hmpc_pcc_comparison_metrics.csv'));
writetable(gain,fullfile('results','hmpc_pcc_comparison_gain.csv'));

category = categorical(stageDisplay,stageDisplay,'Ordinal',true);
fig = figure('Visible','off','Color','w','Position',[45 35 1250 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
bar(category,hmpcMatrix.'); ylim([0 1]); grid on; box off;
ylabel('HMPC'); title('Held-out multi-channel phase concentration');
legend(condition,'Interpreter','none','Location','best');
nexttile;
gainMatrix = [speedGain,pathGain,totalGain];
gainCategory = categorical({'Speed','Path','Total'}, ...
    {'Speed','Path','Total'},'Ordinal',true);
bar(gainCategory,gainMatrix.'); grid on; box off;
yline(0,'k:','HandleVisibility','off');
ylabel('Normalized HMPC gain'); title('Separated phase-alignment gains');
legend(condition,'Interpreter','none','Location','best');
nexttile;
mapRaw = result{2}.stage(1).hmpcMap;
imagesc(linspace(0,360,size(mapRaw,1)),1:size(mapRaw,2),mapRaw.');
axis xy; clim([0 1]); colorbar; colormap(turbo);
xlabel('Mechanical repeat phase (deg)'); ylabel('Mesh harmonic index');
title('PF50 nominal-time phase concentration');
nexttile;
mapFull = result{2}.stage(3).hmpcMap;
imagesc(linspace(0,360,size(mapFull,1)),1:size(mapFull,2),mapFull.');
axis xy; clim([0 1]); colorbar;
xlabel('Mechanical repeat phase (deg)'); ylabel('Mesh harmonic index');
title('PF50 common + cross-fitted path');
exportgraphics(fig,fullfile('results','hmpc_comparison.png'), ...
    'Resolution',190);
close(fig);

fig = figure('Visible','off','Color','w','Position',[55 45 1250 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
bar(category,pccMatrix.'); grid on; box off;
ylabel('PCC'); title('Shi et al. waveform consistency metric');
legend(condition,'Interpreter','none','Location','best');
nexttile;
yyaxis left;
bar(category,reshape(meanCycleStd,nStage,nCondition).');
ylabel('Mean cycle STD');
yyaxis right;
plot(category,reshape(cycleCorrelation,nStage,nCondition).', ...
    'o-','LineWidth',1.1);
ylabel('Cycle correlation'); ylim([-1 1]); grid on; box off;
title('Waveform repeatability components');
nexttile;
rawCycle = result{2}.stage(1).waveformMatrix;
imagesc(linspace(0,360,size(rawCycle,2)), ...
    1:size(rawCycle,1),rawCycle); axis xy; colorbar;
xlabel('Mechanical repeat phase (deg)'); ylabel('Repeat index');
title('PF50 nominal-time cycles');
nexttile;
fullCycle = result{2}.stage(3).waveformMatrix;
imagesc(linspace(0,360,size(fullCycle,2)), ...
    1:size(fullCycle,1),fullCycle); axis xy; colorbar;
xlabel('Mechanical repeat phase (deg)'); ylabel('Repeat index');
title('PF50 common + cross-fitted path cycles');
exportgraphics(fig,fullfile('results','pcc_comparison.png'), ...
    'Resolution',190);
close(fig);

comparison.sourceFile = file;
comparison.condition = condition;
comparison.stage = stage;
comparison.targetFs = targetFs;
comparison.segmentDurationSec = segmentDurationSec;
comparison.repeatCarrierTurns = repeatCarrierTurns;
comparison.samplesPerRepeat = samplesPerRepeat;
comparison.record = result;
comparison.metrics = metrics;
comparison.gain = gain;
comparison.hmpcDefinition = ['Weighted mean resultant length of unit ', ...
    'analytic phasors at the same mechanical-repeat phase, pooled across ', ...
    'held-out repeat cycles, selected sensors, and mesh harmonics.'];
comparison.pccDefinition = ['L2 norm of the cycle-average waveform ', ...
    'divided by the L2 norm of pointwise inter-cycle standard deviation.'];
comparison.caution = ['Only about three 31-carrier-turn repeat cycles ', ...
    'are available in a 50 s record; this run establishes the pipeline ', ...
    'but more repeats are required for final inferential statistics.'];
save(fullfile('results','hmpc_pcc_comparison.mat'), ...
    'comparison','-v7.3');

fprintf('HMPC/PCC comparison complete.\n');
for ic = 1:nCondition
    fprintf('  %s: mesh %.4f Hz, sensors %s, path %s\n', ...
        condition(ic),result{ic}.blindMeshFrequencyHz, ...
        strjoin(string(result{ic}.selectedSensorIndex),'+'), ...
        result{ic}.selectedCandidate);
    fprintf('    cross-fitted amplitude error %.3e\n', ...
        result{ic}.crossFitAmplitudePreservationError);
end
disp(metrics);
disp(gain);
end

function result = processRecord(file,targetFs,startSec,durationSec, ...
    repeatCarrierTurns,samplesPerRepeat,zr)
[signal,sourceFs] = loadSegment(file,startSec,durationSec);
signal = signal-mean(signal,1);
signal = resample(signal,targetFs,sourceFs);
fs = targetFs;
[blindMeshFrequencyHz,spectrum] = blindMeshPeak(signal,fs,[145 190]);
harmonicOrder = 1:4;
centerFrequencyHz = harmonicOrder*blindMeshFrequencyHz;
halfBandwidthHz = 30+4*(harmonicOrder-1);
extracted = pg_extract_analytic_harmonics(signal,fs, ...
    centerFrequencyHz,halfBandwidthHz);
edge = round(0.50*fs);
use = (1+edge):(size(signal,1)-edge);
component = extracted.component(use,:,:);
t = (0:numel(use)-1)'/fs;

options.nIteration = 0;
options.spatialMultiple = 3;
options.pathOrder = 4;
% The paper method uses the physically admissible three-planet orders
% ..., -6, -3, 0, 3, 6, ... .  The unconstrained -12:12 model is retained
% only as an explicit ablation in run_real_phase_ablation.
options.spatialOrderOverride = [];
options.ridge = 2e-3;
options.phaseSmoothSec = 0.012;
options.frequencySmoothSec = 0.10;
common = pg_joint_tacholess_phase(component,fs,harmonicOrder, ...
    blindMeshFrequencyHz,zr,options);

safeOptions.cvStride = 5;
safeOptions.ridge = options.ridge;
safeOptions.minCvGain = 0.015;
safeOptions.maxFoldLoss = 0.010;
safeOptions.stabilityPenalty = 0.25;
safeOptions.returnCrossFittedComponent = true;
safe = pg_select_safe_path_correction(component,common.carrierPhase, ...
    common.spatialOrder,common.pathPenalty,safeOptions);

sensorIndex = safe.selectedSensorIndex;
rawSelected = component(:,sensorIndex,:);
crossFitted = safe.crossFittedCorrectedSelectedComponent;
crossFitAmplitudePreservationError = max(abs(abs(crossFitted(:))- ...
    abs(rawSelected(:))))/max(abs(rawSelected(:)));
nominalMeshPhase = 2*pi*blindMeshFrequencyHz*t;
nominalRepeatPhase = nominalMeshPhase/zr/repeatCarrierTurns;
blindRepeatPhase = common.carrierPhase/repeatCarrierTurns;

stage(1) = evaluateStage(rawSelected,nominalRepeatPhase, ...
    samplesPerRepeat);
stage(2) = evaluateStage(rawSelected,blindRepeatPhase, ...
    samplesPerRepeat);
stage(3) = evaluateStage(crossFitted,blindRepeatPhase, ...
    samplesPerRepeat);

result.sourceFile = file;
result.sourceFs = sourceFs;
result.fs = fs;
result.blindMeshFrequencyHz = blindMeshFrequencyHz;
result.blindSpectrumFrequency = spectrum.frequency;
result.blindSpectrumAmplitude = spectrum.amplitude;
result.harmonicOrder = harmonicOrder;
result.selectedCandidate = safe.selectedCandidate;
result.selectedSensorIndex = sensorIndex;
result.pathCvTable = safe.cvTable;
result.crossFitAmplitudePreservationError = ...
    crossFitAmplitudePreservationError;
result.mechanicalRepeatSec = repeatCarrierTurns*zr/ ...
    blindMeshFrequencyHz;
result.stage = stage;
end

function stage = evaluateStage(component,repeatPhase,samplesPerRepeat)
[nSample,nSensor,nHarmonic] = size(component); %#ok<ASGLU>
waveform = real(sum(mean(component,2),3));
wave = pg_phase_synchronous_matrix(waveform,repeatPhase, ...
    samplesPerRepeat);
waveformMatrix = wave.matrix(:,:,1);
meanWaveform = mean(waveformMatrix,1);
stdWaveform = std(waveformMatrix,1,1);
pcc = norm(meanWaveform)/max(norm(stdWaveform),eps);

hmpcMap = zeros(samplesPerRepeat,nHarmonic);
mapWeight = zeros(samplesPerRepeat,nHarmonic);
for ih = 1:nHarmonic
    phaseMatrix = pg_phase_synchronous_matrix(component(:,:,ih), ...
        repeatPhase,samplesPerRepeat);
    matrix = phaseMatrix.matrix;
    amplitude = abs(matrix);
    weight = zeros(size(amplitude));
    for is = 1:nSensor
        scale = median(amplitude(:,:,is),'all');
        normalized = amplitude(:,:,is)/max(scale,eps);
        weight(:,:,is) = normalized.^2./(1+normalized.^2);
    end
    unit = matrix./max(amplitude,eps);
    numerator = squeeze(sum(sum(weight.*unit,1),3));
    denominator = squeeze(sum(sum(weight,1),3));
    hmpcMap(:,ih) = abs(numerator./max(denominator,eps));
    mapWeight(:,ih) = denominator;
end
hmpc = sum(hmpcMap.*mapWeight,'all')/max(sum(mapWeight,'all'),eps);

stage.hmpc = hmpc;
stage.pcc = pcc;
stage.cycleCorrelation = meanOffDiagonalCorrelation(waveformMatrix);
stage.meanCycleStd = mean(stdWaveform);
stage.cycleCount = size(waveformMatrix,1);
stage.hmpcMap = hmpcMap;
stage.waveformMatrix = waveformMatrix;
end

function gain = normalizedGain(before,after)
gain = (after-before)./max(1-before,eps);
end

function [signal,fs] = loadSegment(file,startSec,durationSec)
variables = whos('-file',file);
name = string({variables.name});
if any(name=="Data")
    dataInfo = variables(name=="Data");
    nSample = dataInfo.size(1);
    fs = 100000;
    loadedRate = load(file,'SampleRate');
    if isfield(loadedRate,'SampleRate')
        fs = double(loadedRate.SampleRate);
    end
    first = round(startSec*fs)+1;
    last = min(nSample,first+round(durationSec*fs)-1);
    try
        source = matfile(file);
        signal = double(source.Data(first:last,1:3));
    catch
        loaded = load(file,'Data');
        signal = double(loaded.Data(first:last,1:3));
    end
elseif any(name=="chanvals")
    fs = 51200;
    dataInfo = variables(name=="chanvals");
    nSample = dataInfo.size(1);
    first = round(startSec*fs)+1;
    last = min(nSample,first+round(durationSec*fs)-1);
    try
        source = matfile(file);
        signal = double(source.chanvals(first:last,1:3));
    catch
        loaded = load(file,'chanvals');
        signal = double(loaded.chanvals(first:last,1:3));
    end
else
    error('run_hmpc_pcc_comparison:Variable', ...
        'Expected Data or chanvals in %s.',file);
end
if size(signal,1)<round(0.90*durationSec*fs)
    error('run_hmpc_pcc_comparison:Duration', ...
        'The record is too short for the requested validation segment.');
end
end

function [frequencyHz,out] = blindMeshPeak(signal,fs,searchBand)
n = size(signal,1);
window = hann(n);
spectrum = 2*abs(fft(signal.*window,[],1))/sum(window);
frequency = (0:n-1)'*fs/n;
combined = sqrt(mean(spectrum.^2,2));
use = frequency>=searchBand(1) & frequency<=searchBand(2);
candidateFrequency = frequency(use);
candidateAmplitude = combined(use);
[~,peak] = max(candidateAmplitude);
frequencyHz = candidateFrequency(peak);
keep = 1:floor(n/2);
out.frequency = frequency(keep);
out.amplitude = spectrum(keep,:);
end

function value = meanOffDiagonalCorrelation(matrix)
if size(matrix,1)<2
    value = nan;
    return;
end
correlation = corrcoef(matrix.');
use = triu(true(size(correlation)),1);
value = mean(correlation(use),'omitnan');
end
