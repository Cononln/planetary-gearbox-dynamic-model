function ablation = run_real_phase_ablation(healthyFile,faultFile)
%RUN_REAL_PHASE_ABLATION Controlled module ablation on measured data.
% MATLAB performs calculations and exports clean source data only.  Paper
% figures are generated separately in Python so all visual outputs share a
% single publication backend.

setup_paths;
rootDir = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || strlength(string(healthyFile))==0
    healthyFile = fullfile(rootDir,'data_local','bl','healthy_record.MAT');
end
if nargin < 2 || strlength(string(faultFile))==0
    faultFile = fullfile(rootDir,'data_local','planet','pf50_record.MAT');
end

sourceFile = {healthyFile;faultFile};
condition = ["BL600_healthy";"PF50_planet_fault"];
method = ["raw_nominal";"common_phase_only"; ...
    "unconstrained_safe_crossfit";"hard3k_all_sensor_fullfit"; ...
    "hard3k_safe_fullfit";"hard3k_safe_crossfit"];
methodDisplay = ["Raw nominal";"Common phase"; ...
    "No spatial constraint";"All sensors, no safeguard"; ...
    "No cross-fitting";"Full method"];
targetFs = 5120;
segmentStartSec = 5;
segmentDurationSec = 50;
repeatCarrierTurns = 31;
samplesPerRepeat = 32768;
ringToothCount = 84;

nCondition = numel(condition);
nMethod = numel(method);
record = cell(nCondition,1);
for ic = 1:nCondition
    record{ic} = processRecord(sourceFile{ic},targetFs, ...
        segmentStartSec,segmentDurationSec,repeatCarrierTurns, ...
        samplesPerRepeat,ringToothCount,method,methodDisplay);
end

nRow = nCondition*nMethod;
conditionColumn = strings(nRow,1);
methodColumn = strings(nRow,1);
methodDisplayColumn = strings(nRow,1);
sensorSubset = strings(nRow,1);
pathTopology = strings(nRow,1);
spatialConstraint = strings(nRow,1);
validationScheme = strings(nRow,1);
blindMeshFrequencyHz = zeros(nRow,1);
mechanicalRepeatSec = zeros(nRow,1);
repeatCycleCount = zeros(nRow,1);
hmpc = zeros(nRow,1);
pcc = zeros(nRow,1);
cycleCorrelation = zeros(nRow,1);
meanCycleStd = zeros(nRow,1);
faultFrequencyHz = zeros(nRow,1);
commonFaultFamilyToBackgroundDb = zeros(nRow,1);
differentialFaultFamilyToBackgroundDb = zeros(nRow,1);
totalEnergyFaultFamilyToBackgroundDb = zeros(nRow,1);
tpsvdRank1EnergyFraction = zeros(nRow,1);
tpsvdRank3EnergyFraction = zeros(nRow,1);
faultCycleRepeatability = zeros(nRow,1);
amplitudePreservationError = zeros(nRow,1);
row = 0;
for ic = 1:nCondition
    for im = 1:nMethod
        row = row+1;
        item = record{ic}.method(im);
        conditionColumn(row) = condition(ic);
        methodColumn(row) = method(im);
        methodDisplayColumn(row) = methodDisplay(im);
        sensorSubset(row) = item.sensorSubset;
        pathTopology(row) = item.pathTopology;
        spatialConstraint(row) = item.spatialConstraint;
        validationScheme(row) = item.validationScheme;
        blindMeshFrequencyHz(row) = record{ic}.blindMeshFrequencyHz;
        mechanicalRepeatSec(row) = record{ic}.mechanicalRepeatSec;
        repeatCycleCount(row) = item.phaseMetric.cycleCount;
        hmpc(row) = item.phaseMetric.hmpc;
        pcc(row) = item.phaseMetric.pcc;
        cycleCorrelation(row) = item.phaseMetric.cycleCorrelation;
        meanCycleStd(row) = item.phaseMetric.meanCycleStd;
        faultFrequencyHz(row) = record{ic}.faultFrequencyHz;
        commonFaultFamilyToBackgroundDb(row) = ...
            item.diagnosis.common.faultFamilyToBackgroundDb;
        differentialFaultFamilyToBackgroundDb(row) = ...
            item.diagnosis.differential.faultFamilyToBackgroundDb;
        totalEnergyFaultFamilyToBackgroundDb(row) = ...
            item.diagnosis.totalEnergy.faultFamilyToBackgroundDb;
        tpsvdRank1EnergyFraction(row) = ...
            item.diagnosis.totalEnergy.tpsvdRank1EnergyFraction;
        tpsvdRank3EnergyFraction(row) = ...
            item.diagnosis.totalEnergy.tpsvdRank3EnergyFraction;
        faultCycleRepeatability(row) = ...
            item.diagnosis.totalEnergy.faultCycleRepeatability;
        amplitudePreservationError(row) = ...
            item.amplitudePreservationError;
    end
end
metrics = table(conditionColumn,methodColumn,methodDisplayColumn, ...
    sensorSubset,pathTopology,spatialConstraint,validationScheme, ...
    blindMeshFrequencyHz,mechanicalRepeatSec,repeatCycleCount,hmpc,pcc, ...
    cycleCorrelation,meanCycleStd,faultFrequencyHz, ...
    commonFaultFamilyToBackgroundDb,differentialFaultFamilyToBackgroundDb, ...
    totalEnergyFaultFamilyToBackgroundDb,tpsvdRank1EnergyFraction, ...
    tpsvdRank3EnergyFraction,faultCycleRepeatability, ...
    amplitudePreservationError, ...
    'VariableNames',{'condition','method','method_display', ...
    'sensor_subset','path_topology','spatial_constraint', ...
    'validation_scheme','blind_mesh_frequency_hz', ...
    'mechanical_repeat_sec','repeat_cycle_count','hmpc','pcc', ...
    'cycle_correlation','mean_cycle_std','planet_fault_frequency_hz', ...
    'common_fault_family_to_background_db', ...
    'differential_fault_family_to_background_db', ...
    'total_energy_fault_family_to_background_db', ...
    'tpsvd_rank1_energy_fraction', ...
    'tpsvd_rank3_energy_fraction','fault_cycle_repeatability', ...
    'amplitude_preservation_error'});

metricMatrix = @(value) reshape(value,nMethod,nCondition).';
hmpcMatrix = metricMatrix(hmpc);
pccMatrix = metricMatrix(pcc);
commonFaultSnrMatrix = metricMatrix(commonFaultFamilyToBackgroundDb);
differentialFaultSnrMatrix = ...
    metricMatrix(differentialFaultFamilyToBackgroundDb);
totalFaultSnrMatrix = metricMatrix(totalEnergyFaultFamilyToBackgroundDb);
rank3Matrix = metricMatrix(tpsvdRank3EnergyFraction);
repeatabilityMatrix = metricMatrix(faultCycleRepeatability);
contrast = table(method,methodDisplay, ...
    (hmpcMatrix(2,:)-hmpcMatrix(1,:)).', ...
    (pccMatrix(2,:)-pccMatrix(1,:)).', ...
    (commonFaultSnrMatrix(2,:)-commonFaultSnrMatrix(1,:)).', ...
    (differentialFaultSnrMatrix(2,:)- ...
    differentialFaultSnrMatrix(1,:)).', ...
    (totalFaultSnrMatrix(2,:)-totalFaultSnrMatrix(1,:)).', ...
    (rank3Matrix(2,:)-rank3Matrix(1,:)).', ...
    (repeatabilityMatrix(2,:)-repeatabilityMatrix(1,:)).', ...
    'VariableNames',{'method','method_display','PF50_minus_BL_hmpc', ...
    'PF50_minus_BL_pcc','PF50_minus_BL_common_fault_snr_db', ...
    'PF50_minus_BL_differential_fault_snr_db', ...
    'PF50_minus_BL_total_energy_fault_snr_db', ...
    'PF50_minus_BL_tpsvd_rank3_fraction', ...
    'PF50_minus_BL_fault_cycle_repeatability'});

if ~isfolder('results')
    mkdir('results');
end
writetable(metrics,fullfile('results','phase_ablation_metrics.csv'));
writetable(contrast,fullfile('results','phase_ablation_contrast.csv'));

spectrumTable = table;
componentName = ["common";"differential";"total_energy"];
for ic = 1:nCondition
    for im = 1:nMethod
        item = record{ic}.method(im);
        diagnosisSet = {item.diagnosis.common; ...
            item.diagnosis.differential;item.diagnosis.totalEnergy};
        for id = 1:numel(diagnosisSet)
            diagnosisItem = diagnosisSet{id};
            keep = diagnosisItem.frequencyHz<=60;
            n = nnz(keep);
            block = table(repmat(condition(ic),n,1), ...
                repmat(method(im),n,1),repmat(methodDisplay(im),n,1), ...
                repmat(componentName(id),n,1), ...
                diagnosisItem.frequencyHz(keep), ...
                diagnosisItem.spectrumAmplitude(keep), ...
                diagnosisItem.spectrumRelativeDb(keep), ...
                'VariableNames',{'condition','method','method_display', ...
                'signal_component','frequency_hz','amplitude', ...
                'relative_to_background_db'});
            spectrumTable = [spectrumTable;block]; %#ok<AGROW>
        end
    end
end
writetable(spectrumTable,fullfile('results', ...
    'phase_ablation_spectra_source.csv'));

mapTable = makeMapSource(record{2},method([1,6]),methodDisplay([1,6]));
cycleTable = makeCycleSource(record{2},method([1,6]), ...
    methodDisplay([1,6]));
writetable(mapTable,fullfile('results','phase_ablation_hmpc_map_source.csv'));
writetable(cycleTable,fullfile('results','phase_ablation_cycle_source.csv'));

ablation.sourceFile = sourceFile;
ablation.condition = condition;
ablation.methodName = method;
ablation.methodDisplay = methodDisplay;
ablation.targetFs = targetFs;
ablation.segmentDurationSec = segmentDurationSec;
ablation.repeatCarrierTurns = repeatCarrierTurns;
ablation.samplesPerRepeat = samplesPerRepeat;
ablation.record = record;
ablation.metrics = metrics;
ablation.contrast = contrast;
ablation.metricDefinitions = struct( ...
    'HMPC',['Amplitude-reliability-weighted resultant length of unit ', ...
    'analytic phasors pooled across held-out repeat cycles, sensors, ', ...
    'and mesh harmonics.'], ...
    'PCC',['L2 norm of the cycle-average waveform divided by the L2 ', ...
    'norm of pointwise inter-cycle standard deviation.'], ...
    'faultFamilySNR',['RMS of the first four theoretical planet-fault ', ...
    'spectral lines divided by the median 2-60 Hz off-line background; ', ...
    'reported separately for aligned common, differential residual, and ', ...
    'energy-preserving total envelopes.']);
ablation.statisticsCaution = ['Each 50 s record contains only three ', ...
    'complete 31-carrier-turn mechanical repeats.  Results are controlled ', ...
    'within-record ablations; no independent-record significance test is ', ...
    'claimed.'];
save(fullfile('results','phase_ablation.mat'),'ablation','-v7.3');

fprintf('Measured-data phase ablation complete.\n');
for ic = 1:nCondition
    fprintf('  %s: mesh %.4f Hz, full-method sensors %s, path %s\n', ...
        condition(ic),record{ic}.blindMeshFrequencyHz, ...
        record{ic}.method(6).sensorSubset, ...
        record{ic}.method(6).pathTopology);
end
disp(metrics);
disp(contrast);
end

function result = processRecord(file,targetFs,startSec,durationSec, ...
    repeatCarrierTurns,samplesPerRepeat,zr,methodName,methodDisplay)
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

commonOptions.nIteration = 0;
commonOptions.spatialMultiple = 3;
commonOptions.pathOrder = 4;
commonOptions.ridge = 2e-3;
commonOptions.phaseSmoothSec = 0.012;
commonOptions.frequencySmoothSec = 0.10;
commonOptions.spatialOrderOverride = [];
commonOptions.offGridPenalty = 1;
common = pg_joint_tacholess_phase(component,fs,harmonicOrder, ...
    blindMeshFrequencyHz,zr,commonOptions);

safeOptions.cvStride = 5;
safeOptions.ridge = commonOptions.ridge;
safeOptions.minCvGain = 0.015;
safeOptions.maxFoldLoss = 0.010;
safeOptions.stabilityPenalty = 0.25;
safeOptions.returnCrossFittedComponent = true;
hardSafe = pg_select_safe_path_correction(component, ...
    common.carrierPhase,common.spatialOrder,common.pathPenalty,safeOptions);

freeOrder = -12:12;
freePenalty = ones(size(freeOrder));
freeSafe = pg_select_safe_path_correction(component, ...
    common.carrierPhase,freeOrder,freePenalty,safeOptions);

selectedIndex = hardSafe.selectedSensorIndex;
selectedRaw = component(:,selectedIndex,:);
freeIndex = freeSafe.selectedSensorIndex;
freeRaw = component(:,freeIndex,:);
nominalMeshPhase = 2*pi*blindMeshFrequencyHz*t;
nominalRepeatPhase = nominalMeshPhase/zr/repeatCarrierTurns;
blindRepeatPhase = common.carrierPhase/repeatCarrierTurns;

componentSet = {selectedRaw,selectedRaw, ...
    freeSafe.crossFittedCorrectedSelectedComponent, ...
    common.correctedComponent,hardSafe.correctedSelectedComponent, ...
    hardSafe.crossFittedCorrectedSelectedComponent};
referenceSet = {selectedRaw,selectedRaw,freeRaw,component, ...
    selectedRaw,selectedRaw};
repeatPhaseSet = {nominalRepeatPhase,blindRepeatPhase,blindRepeatPhase, ...
    blindRepeatPhase,blindRepeatPhase,blindRepeatPhase};
meshPhaseSet = {nominalMeshPhase,common.meshPhase,common.meshPhase, ...
    common.meshPhase,common.meshPhase,common.meshPhase};
sensorIndexSet = {selectedIndex,selectedIndex,freeIndex,1:3, ...
    selectedIndex,selectedIndex};
pathTopology = ["none";"none";freeSafe.selectedCandidate; ...
    "complete_graph";hardSafe.selectedCandidate;hardSafe.selectedCandidate];
spatialConstraint = ["none";"none";"orders_-12_to_12"; ...
    "three_planet_3k";"three_planet_3k";"three_planet_3k"];
validationScheme = ["none";"none";"cross_fitted";"full_fit"; ...
    "full_fit";"cross_fitted"];

nMethod = numel(methodName);
method = repmat(struct, nMethod,1);
faultFrequencyHz = 2*blindMeshFrequencyHz/31;
for im = 1:nMethod
    method(im).name = methodName(im);
    method(im).display = methodDisplay(im);
    method(im).sensorSubset = strjoin(string(sensorIndexSet{im}),'+');
    method(im).pathTopology = pathTopology(im);
    method(im).spatialConstraint = spatialConstraint(im);
    method(im).validationScheme = validationScheme(im);
    method(im).amplitudePreservationError = amplitudeError( ...
        componentSet{im},referenceSet{im});
    method(im).phaseMetric = evaluatePhaseMetric(componentSet{im}, ...
        repeatPhaseSet{im},samplesPerRepeat);
    [commonEnvelope,differentialEnvelope,totalEnergyEnvelope] = ...
        separatedEnvelopes(componentSet{im},referenceSet{im});
    method(im).diagnosis.common = evaluateDiagnosis(commonEnvelope, ...
        meshPhaseSet{im},fs,faultFrequencyHz);
    method(im).diagnosis.differential = evaluateDiagnosis( ...
        differentialEnvelope,meshPhaseSet{im},fs,faultFrequencyHz);
    method(im).diagnosis.totalEnergy = evaluateDiagnosis( ...
        totalEnergyEnvelope,meshPhaseSet{im},fs,faultFrequencyHz);
end

result.sourceFile = file;
result.sourceFs = sourceFs;
result.fs = fs;
result.blindMeshFrequencyHz = blindMeshFrequencyHz;
result.blindSpectrumFrequency = spectrum.frequency;
result.blindSpectrumAmplitude = spectrum.amplitude;
result.harmonicOrder = harmonicOrder;
result.mechanicalRepeatSec = repeatCarrierTurns*zr/ ...
    blindMeshFrequencyHz;
result.faultFrequencyHz = faultFrequencyHz;
result.hardSafeCvTable = hardSafe.cvTable;
result.freeSafeCvTable = freeSafe.cvTable;
result.method = method;
end

function phaseMetric = evaluatePhaseMetric(component,repeatPhase, ...
    samplesPerRepeat)
[~,nSensor,nHarmonic] = size(component);
waveform = real(sum(mean(component,2),3));
wave = pg_phase_synchronous_matrix(waveform,repeatPhase,samplesPerRepeat);
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

phaseMetric.hmpc = hmpc;
phaseMetric.pcc = pcc;
phaseMetric.cycleCorrelation = meanOffDiagonalCorrelation(waveformMatrix);
phaseMetric.meanCycleStd = mean(stdWaveform);
phaseMetric.cycleCount = size(waveformMatrix,1);
phaseMetric.hmpcMap = hmpcMap;
phaseMetric.waveformMatrix = waveformMatrix;
end

function [commonEnvelope,differentialEnvelope,totalEnergyEnvelope] = ...
    separatedEnvelopes(component,reference)
[nSample,~,nHarmonic] = size(component);
commonEnvelope = zeros(nSample,1);
differentialEnvelope = zeros(nSample,1);
totalEnergyEnvelope = zeros(nSample,1);
for ih = 1:nHarmonic
    scale = sqrt(mean(abs(reference(:,:,ih)).^2,'all'));
    common = mean(component(:,:,ih),2);
    residual = component(:,:,ih)-common;
    residualRms = sqrt(mean(abs(residual).^2,2));
    commonEnvelope = commonEnvelope+abs(common)/max(scale,eps);
    differentialEnvelope = differentialEnvelope+ ...
        residualRms/max(scale,eps);
    totalEnergyEnvelope = totalEnergyEnvelope+ ...
        sqrt(abs(common).^2+residualRms.^2)/ ...
        max(scale,eps);
end
commonEnvelope = detrend(commonEnvelope/nHarmonic);
differentialEnvelope = detrend(differentialEnvelope/nHarmonic);
totalEnergyEnvelope = detrend(totalEnergyEnvelope/nHarmonic);
end

function diagnosis = evaluateDiagnosis(envelope,meshPhase,fs, ...
    faultFrequencyHz)
[frequencyHz,spectrumAmplitude] = oneSidedSpectrum(envelope,fs);
lineAmplitude = zeros(4,1);
background = frequencyHz>=2 & frequencyHz<=60;
for ih = 1:4
    target = ih*faultFrequencyHz;
    useLine = abs(frequencyHz-target)<=0.35;
    lineAmplitude(ih) = max(spectrumAmplitude(useLine));
    background = background & abs(frequencyHz-target)>0.9;
end
backgroundLevel = median(spectrumAmplitude(background));
lineRms = sqrt(mean(lineAmplitude.^2));
faultFamilyToBackgroundDb = 20*log10(lineRms/max(backgroundLevel,eps));
spectrumRelativeDb = 20*log10(max(spectrumAmplitude,eps)/ ...
    max(backgroundLevel,eps));

faultPhase = 2*meshPhase/31;
phaseMatrix = pg_phase_synchronous_matrix(envelope,faultPhase,512);
matrix = phaseMatrix.matrix(:,:,1);
matrix = matrix-mean(matrix,2);
[~,s,~] = svd(matrix,'econ');
singularValue = diag(s);
rank1 = singularValue(1)^2/sum(singularValue.^2);
rank3 = sum(singularValue(1:min(3,numel(singularValue))).^2)/ ...
    sum(singularValue.^2);

diagnosis.frequencyHz = frequencyHz;
diagnosis.spectrumAmplitude = spectrumAmplitude;
diagnosis.spectrumRelativeDb = spectrumRelativeDb;
diagnosis.backgroundLevel = backgroundLevel;
diagnosis.lineAmplitude = lineAmplitude;
diagnosis.faultFamilyToBackgroundDb = faultFamilyToBackgroundDb;
diagnosis.tpsvdRank1EnergyFraction = rank1;
diagnosis.tpsvdRank3EnergyFraction = rank3;
diagnosis.faultCycleRepeatability = ...
    meanOffDiagonalCorrelation(matrix);
diagnosis.faultCycleCount = size(matrix,1);
diagnosis.singularValue = singularValue;
end

function mapTable = makeMapSource(faultRecord,methodName,methodDisplay)
mapTable = table;
methodIndex = [1,6];
for k = 1:numel(methodIndex)
    im = methodIndex(k);
    map = faultRecord.method(im).phaseMetric.hmpcMap;
    sampleIndex = unique(round(linspace(1,size(map,1),720)));
    phaseDeg = (sampleIndex-1)/(size(map,1)-1)*360;
    for ih = 1:size(map,2)
        n = numel(sampleIndex);
        block = table(repmat(methodName(k),n,1), ...
            repmat(methodDisplay(k),n,1),phaseDeg(:), ...
            repmat(ih,n,1),map(sampleIndex,ih), ...
            'VariableNames',{'method','method_display', ...
            'mechanical_repeat_phase_deg','mesh_harmonic', ...
            'phase_concentration'});
        mapTable = [mapTable;block]; %#ok<AGROW>
    end
end
end

function cycleTable = makeCycleSource(faultRecord,methodName,methodDisplay)
cycleTable = table;
methodIndex = [1,6];
for k = 1:numel(methodIndex)
    im = methodIndex(k);
    matrix = faultRecord.method(im).phaseMetric.waveformMatrix;
    matrix = matrix/max(sqrt(mean(matrix.^2,'all')),eps);
    sampleIndex = unique(round(linspace(1,size(matrix,2),720)));
    phaseDeg = (sampleIndex-1)/(size(matrix,2)-1)*360;
    for ir = 1:size(matrix,1)
        n = numel(sampleIndex);
        block = table(repmat(methodName(k),n,1), ...
            repmat(methodDisplay(k),n,1),repmat(ir,n,1), ...
            phaseDeg(:),matrix(ir,sampleIndex).', ...
            'VariableNames',{'method','method_display','repeat_index', ...
            'mechanical_repeat_phase_deg','normalized_amplitude'});
        cycleTable = [cycleTable;block]; %#ok<AGROW>
    end
end
end

function value = amplitudeError(component,reference)
value = max(abs(abs(component(:))-abs(reference(:))))/ ...
    max(abs(reference(:)));
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
    error('run_real_phase_ablation:Variable', ...
        'Expected Data or chanvals in %s.',file);
end
if size(signal,1)<round(0.90*durationSec*fs)
    error('run_real_phase_ablation:Duration', ...
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

function [frequency,amplitude] = oneSidedSpectrum(signal,fs)
n = numel(signal);
window = hann(n);
fullAmplitude = 2*abs(fft(signal.*window))/sum(window);
nHalf = floor(n/2);
amplitude = fullAmplitude(1:nHalf);
frequency = (0:nHalf-1)'*fs/n;
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
