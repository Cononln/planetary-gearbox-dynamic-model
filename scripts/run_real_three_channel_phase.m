function report = run_real_three_channel_phase(matFile)
%RUN_REAL_THREE_CHANNEL_PHASE Blind three-channel phase decomposition.
% The measured vibration channels are located at 0, 120, and 240 degrees.
% For legacy chanvals files, columns 4/5 are quadrature carrier-encoder
% channels.  They remain untouched until every vibration-only estimate and
% regularization choice has been frozen.

setup_paths;
rootDir = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || strlength(string(matFile))==0
    matFile = fullfile(rootDir,'data_local','planet','P050S600L00.tdms.mat');
end

zs = 21;
zr = 84;
carrierEncoderPpr = 1024;
sensorAngleDeg = [0,120,240];
segmentDurationSec = 12;

[rawSignal,encoderA,fs,segmentIndex,sourceInfo] = ...
    loadThreeChannelSegment(matFile,segmentDurationSec);
rawSignal = rawSignal-mean(rawSignal,1);
tSegment = (0:size(rawSignal,1)-1)'/fs;

expectedMeshFrequencyHz = zs*zr/(zs+zr)*(sourceInfo.nominalRpm/60);
searchBandHz = expectedMeshFrequencyHz+[ -0.12,0.12]* ...
    expectedMeshFrequencyHz;
[blindMeshFrequencyHz,blindSpectrum] = blindMeshPeak( ...
    rawSignal,fs,searchBandHz);

[harmonicOrder,harmonicProminenceTable] = selectReliableHarmonics( ...
    blindSpectrum,blindMeshFrequencyHz,fs,4);
centerFrequencyHz = harmonicOrder*blindMeshFrequencyHz;
halfBandwidthHz = 30+4*(harmonicOrder-1);
extracted = pg_extract_analytic_harmonics(rawSignal,fs, ...
    centerFrequencyHz,halfBandwidthHz);
edgeSec = 0.50;
use = tSegment>=edgeSec & tSegment<=tSegment(end)-edgeSec;
t = tSegment(use)-tSegment(find(use,1));
component = extracted.component(use,:,:);
signalUse = rawSignal(use,:);

baseOptions.nIteration = 3;
baseOptions.spatialMultiple = 3;
baseOptions.pathOrder = 4;
baseOptions.ridge = 2e-3;
baseOptions.relaxation = 0.35;
baseOptions.commonAnchorWeight = 0.80;
baseOptions.phaseSmoothSec = 0.012;
baseOptions.frequencySmoothSec = 0.10;
baseOptions.removeCommonSpatialOrder = false;

% All objects above and all choices below use vibration only.
noPathPhase = blindCommonPhase(component,fs,harmonicOrder, ...
    blindMeshFrequencyHz,baseOptions.phaseSmoothSec);
hardThree = pg_joint_tacholess_phase(component,fs,harmonicOrder, ...
    blindMeshFrequencyHz,zr,baseOptions);

freeOptions = baseOptions;
freeOptions.spatialOrderOverride = -12:12;
freeOptions.offGridPenalty = 1;
freeThree = pg_joint_tacholess_phase(component,fs,harmonicOrder, ...
    blindMeshFrequencyHz,zr,freeOptions);

penaltyCandidate = [1,3,10,30,100,300];
[selectedOffGridPenalty,penaltyCvTable] = selectOffGridPenaltyMulti( ...
    component,hardThree.carrierPhase,12,3,baseOptions.ridge, ...
    penaltyCandidate,freeThree.pathGraph.edgeFitCoherence);
adaptiveOptions = freeOptions;
adaptiveOptions.offGridPenalty = selectedOffGridPenalty;
adaptiveOptions.minEdgeFitCoherence = 0.50;
adaptiveOptions.nIteration = 0;
[selectedCommonSensorMask,commonSensorCvTable] = ...
    selectCommonSensorSubset(component,fs,harmonicOrder, ...
    blindMeshFrequencyHz,baseOptions.phaseSmoothSec);
adaptiveOptions.commonSensorHarmonicWeight = repmat( ...
    selectedCommonSensorMask(:),1,size(component,3));
adaptiveOptions.pathSensorHarmonicWeight = sensorReliabilityFromGraph( ...
    freeThree.pathGraph,size(component,2),size(component,3));
adaptiveThree = pg_joint_tacholess_phase(component,fs,harmonicOrder, ...
    blindMeshFrequencyHz,zr,adaptiveOptions);
adaptiveTwoOptions = adaptiveOptions;
adaptiveTwoOptions.commonSensorHarmonicWeight = ...
    ones(2,size(component,3));
adaptiveTwoOptions.pathSensorHarmonicWeight = ...
    adaptiveOptions.pathSensorHarmonicWeight(1:2,:);
adaptiveTwo = pg_joint_tacholess_phase(component(:,1:2,:),fs, ...
    harmonicOrder,blindMeshFrequencyHz,zr,adaptiveTwoOptions);

% Path topology is selected without an encoder.  Odd/even carrier
% revolutions are held out in turn; an unhelpful three-channel correction
% therefore falls back to the best reliable pair or to no correction.
safeOptions.cvStride = 20;
safeOptions.ridge = baseOptions.ridge;
safeOptions.minCvGain = 0.015;
safeOptions.maxFoldLoss = 0.010;
safeOptions.stabilityPenalty = 0.25;
safePath = pg_select_safe_path_correction(component, ...
    adaptiveThree.carrierPhase,adaptiveThree.spatialOrder, ...
    adaptiveThree.pathPenalty,safeOptions);

%% Encoder truth: intentionally decoded only after blind estimates freeze.
truth = struct;
hasEncoderTruth = ~isempty(encoderA);
if hasEncoderTruth
    encoderSample = encoderA(1:100:end);
    threshold = 0.5*(percentile(encoderSample,1)+ ...
        percentile(encoderSample,99));
    minimumSeparation = max(2,round(fs/(sourceInfo.nominalRpm/60/5* ...
        carrierEncoderPpr)*0.25));
    carrierEdges = risingEdges(encoderA,threshold,minimumSeparation);
    carrierEdgePhase = (0:numel(carrierEdges)-1)'* ...
        (2*pi/carrierEncoderPpr);
    carrierPhaseSegment = interp1(carrierEdges,carrierEdgePhase, ...
        segmentIndex,'linear','extrap');
    truth.carrierPhase = carrierPhaseSegment(use);
    truth.meshPhase = zr*truth.carrierPhase;
    truth.meshFrequencyHz = instantaneousFrequency( ...
        truth.meshPhase,fs,0.10);
    truth.carrierEdges = carrierEdges;
    truth.definition = ['Carrier encoder A-phase, 1024 rising edges per ', ...
        'carrier revolution; withheld until scoring.'];
    oraclePath = fitKnownCarrierPath(component,truth.carrierPhase, ...
        -12:12,baseOptions.ridge);
else
    oraclePath = [];
end

method = {'no_path_three_channel';'reliable_two_channel'; ...
    'hard_3k_three_channel';'free_three_channel'; ...
    'reliability_graph_three_channel'; ...
    'cross_validated_safe_path'};
nSensorUsed = [3;2;3;3;3;numel(safePath.selectedSensorIndex)];
sensorIndexSet = {1:3,1:2,1:3,1:3,1:3, ...
    safePath.selectedSensorIndex};
phaseSet = {noPathPhase,adaptiveTwo.meshPhase,hardThree.meshPhase, ...
    freeThree.meshPhase,adaptiveThree.meshPhase,adaptiveThree.meshPhase};
componentSet = {component,adaptiveTwo.correctedComponent, ...
    hardThree.correctedComponent,freeThree.correctedComponent, ...
    adaptiveThree.correctedComponent, ...
    safePath.correctedSelectedComponent};
pathSet = {[],adaptiveTwo.pathPhasor,hardThree.pathPhasor, ...
    freeThree.pathPhasor,adaptiveThree.pathPhasor, ...
    safePath.selectedPathPhasor};
graphSet = {[],adaptiveTwo.pathGraph,hardThree.pathGraph, ...
    freeThree.pathGraph,adaptiveThree.pathGraph,safePath.pathGraph};
if safePath.selectedCandidate=="no_path"
    pathSet{end} = [];
end

nMethod = numel(method);
meshPhaseRmseDeg = nan(nMethod,1);
meshPhaseP95Deg = nan(nMethod,1);
carrierSpeedRmseRpm = nan(nMethod,1);
pathPhaseVsEncoderRmseDeg = nan(nMethod,1);
meanPairResidualDeg = zeros(nMethod,1);
meanPairCoherence = zeros(nMethod,1);
activePairCoherence = zeros(nMethod,1);
activeSensorFraction = zeros(nMethod,1);
amplitudePreservationError = zeros(nMethod,1);
cycleClosureRmseDeg = nan(nMethod,1);
for im = 1:nMethod
    sensorIndex = sensorIndexSet{im};
    original = component(:,sensorIndex,:);
    estimate = componentSet{im};
    meanPairResidualDeg(im) = pairResidual(estimate);
    meanPairCoherence(im) = pairCoherence(estimate);
    activePairCoherence(im) = pairCoherence(estimate);
    amplitudePreservationError(im) = max(abs(abs(estimate(:))- ...
        abs(original(:))))/max(abs(original(:)));
    if ~isempty(graphSet{im})
        cycleClosureRmseDeg(im) = ...
            graphSet{im}.meanCycleClosureRmseDeg;
        activePairCoherence(im) = graphPairCoherence( ...
            estimate,graphSet{im});
        activeSensorFraction(im) = mean( ...
            graphSet{im}.activeSensorMask,'all');
    else
        activeSensorFraction(im) = 1;
    end
    if hasEncoderTruth
        [meshPhaseRmseDeg(im),meshPhaseP95Deg(im), ...
            carrierSpeedRmseRpm(im)] = phaseScore( ...
            phaseSet{im},truth,fs,zr);
        if ~isempty(pathSet{im})
            pathPhaseVsEncoderRmseDeg(im) = pathScore( ...
                pathSet{im},oraclePath(:,sensorIndex,:),original, ...
                graphSet{im});
        end
    end
end

sensorSubset = strings(nMethod,1);
for im = 1:nMethod
    sensorSubset(im) = strjoin(string(sensorIndexSet{im}),'+');
end
metrics = table(method,nSensorUsed,sensorSubset, ...
    meshPhaseRmseDeg,meshPhaseP95Deg, ...
    carrierSpeedRmseRpm,pathPhaseVsEncoderRmseDeg, ...
    meanPairResidualDeg,meanPairCoherence,activePairCoherence, ...
    activeSensorFraction,cycleClosureRmseDeg, ...
    amplitudePreservationError, ...
    'VariableNames',{'method','sensor_count','sensor_subset', ...
    'mesh_phase_rmse_deg', ...
    'mesh_phase_p95_deg','carrier_speed_rmse_rpm', ...
    'path_phase_vs_encoder_rmse_deg','mean_pair_residual_deg', ...
    'mean_pair_coherence','active_pair_coherence', ...
    'active_sensor_fraction','fitted_cycle_closure_rmse_deg', ...
    'amplitude_preservation_error'});

tag = regexprep(sourceInfo.baseName,'[^A-Za-z0-9]+','_');
prefix = fullfile('results',['three_channel_',tag]);
writetable(metrics,[prefix,'_metrics.csv']);
writetable(penaltyCvTable,[prefix,'_penalty_cv.csv']);
writetable(harmonicProminenceTable,[prefix,'_harmonic_selection.csv']);
writetable(commonSensorCvTable,[prefix,'_common_sensor_cv.csv']);
writetable(safePath.cvTable,[prefix,'_path_safeguard_cv.csv']);

fig = figure('Visible','off','Color','w','Position',[50 40 1240 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
plot(blindSpectrum.frequency,blindSpectrum.amplitude,'LineWidth',0.8);
xlim(searchBandHz); grid on; box off;
xlabel('Frequency (Hz)'); ylabel('Amplitude');
title(sprintf('Three-channel blind mesh peak: %.3f Hz', ...
    blindMeshFrequencyHz));
legend('0 deg','120 deg','240 deg','Location','best');

nexttile;
if hasEncoderTruth
    plot(t,movingAverage(truth.meshFrequencyHz,round(0.50*fs)), ...
        'k','LineWidth',1.2); hold on;
    plot(t,movingAverage(adaptiveTwo.instantaneousMeshFrequencyHz, ...
        round(0.50*fs)),'Color',[0.25 0.55 0.85],'LineWidth',0.8);
    plot(t,movingAverage(adaptiveThree.instantaneousMeshFrequencyHz, ...
        round(0.50*fs)),'Color',[0.85 0.20 0.10],'LineWidth',0.9);
    legend('Hidden encoder','Two channel','Three channel', ...
        'Location','best');
else
    plot(t,adaptiveTwo.instantaneousMeshFrequencyHz,'LineWidth',0.8); hold on;
    plot(t,adaptiveThree.instantaneousMeshFrequencyHz,'LineWidth',0.8);
    legend('Two channel','Three channel','Location','best');
end
xlim([t(1),t(end)]); grid on; box off;
xlabel('Time (s)'); ylabel('Mesh frequency (Hz)');
title('Encoder-hidden common phase estimate');

nexttile;
methodCategory = categorical(method,method,'Ordinal',true);
yyaxis left;
bar(methodCategory,meanPairResidualDeg,'FaceColor',[0.25 0.55 0.80]);
ylabel('Mean pair residual (deg)');
yyaxis right;
plot(methodCategory,meanPairCoherence,'ro-', ...
    'LineWidth',1.2,'MarkerFaceColor','r');
ylabel('Mean pair coherence'); ylim([0 1.05]);
grid on; box off; title('All-pair differential alignment');

nexttile;
if hasEncoderTruth
    yyaxis left;
    bar(methodCategory,meshPhaseRmseDeg,'FaceColor',[0.35 0.65 0.45]);
    ylabel('Mesh-phase RMSE (deg)');
    yyaxis right;
    plot(methodCategory,pathPhaseVsEncoderRmseDeg,'mo-', ...
        'LineWidth',1.2,'MarkerFaceColor','m');
    ylabel('Path-phase RMSE (deg)');
    title('Hidden-encoder scoring');
else
    bar(methodCategory,cycleClosureRmseDeg, ...
        'FaceColor',[0.55 0.40 0.75]);
    ylabel('Fitted closure RMSE (deg)');
    title('Three-channel fitted-path consistency');
end
grid on; box off;
exportgraphics(fig,[prefix,'_validation.png'],'Resolution',190);
close(fig);

if hasEncoderTruth
    fig = figure('Visible','off','Color','w','Position',[70 60 1240 430]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    pair = [1,2;1,3;2,3];
    ihPlot = numel(harmonicOrder);
    for ip = 1:3
        nexttile;
        oracleDifference = oraclePath(:,pair(ip,2),ihPlot).* ...
            conj(oraclePath(:,pair(ip,1),ihPlot));
        blindDifference = adaptiveThree.pathPhasor(:,pair(ip,2),ihPlot).* ...
            conj(adaptiveThree.pathPhasor(:,pair(ip,1),ihPlot));
        [angleDeg,oracleMap] = binnedPhaseMap( ...
            truth.carrierPhase,oracleDifference,180);
        [~,blindMap] = binnedPhaseMap( ...
            truth.carrierPhase,blindDifference,180);
        plot(angleDeg,oracleMap,'k','LineWidth',1.1); hold on;
        plot(angleDeg,blindMap,'Color',[0.85 0.20 0.10], ...
            'LineWidth',0.8);
        xlim([0 360]); ylim([-180 180]); grid on; box off;
        xlabel('Carrier angle (deg)'); ylabel('Path phase (deg)');
        title(sprintf('%d deg - %d deg, h=%d', ...
            sensorAngleDeg(pair(ip,2)),sensorAngleDeg(pair(ip,1)), ...
            harmonicOrder(ihPlot)));
    end
    legend('Encoder-angle oracle','Blind graph estimate', ...
        'Location','best');
    exportgraphics(fig,[prefix,'_path_maps.png'],'Resolution',190);
    close(fig);
end

report.sourceFile = matFile;
report.sourceInfo = sourceInfo;
report.fs = fs;
report.sensorAngleDeg = sensorAngleDeg;
report.segmentIndex = segmentIndex;
report.t = t;
report.rawSignal = signalUse;
report.blindMeshFrequencyHz = blindMeshFrequencyHz;
report.harmonicOrder = harmonicOrder;
report.harmonicProminenceTable = harmonicProminenceTable;
report.selectedOffGridPenalty = selectedOffGridPenalty;
report.penaltyCvTable = penaltyCvTable;
report.selectedCommonSensorMask = selectedCommonSensorMask;
report.commonSensorCvTable = commonSensorCvTable;
report.hardThree = compactPhaseResult(hardThree);
report.freeThree = compactPhaseResult(freeThree);
report.adaptiveTwo = compactPhaseResult(adaptiveTwo);
report.adaptiveThree = compactPhaseResult(adaptiveThree);
report.adaptiveThree.correctedSignal = adaptiveThree.correctedSignal;
report.adaptiveThree.pathPhasorDecimated = ...
    adaptiveThree.pathPhasor(1:20:end,:,:);
report.safePath.selectedCandidate = safePath.selectedCandidate;
report.safePath.selectedSensorMask = safePath.selectedSensorMask;
report.safePath.selectedSensorIndex = safePath.selectedSensorIndex;
report.safePath.selectedPairMask = safePath.selectedPairMask;
report.safePath.cvTable = safePath.cvTable;
report.safePath.fullRawCoherence = safePath.fullRawCoherence;
report.safePath.fullCorrectedCoherence = ...
    safePath.fullCorrectedCoherence;
report.safePath.fullCoherenceGain = safePath.fullCoherenceGain;
report.safePath.amplitudePreservationError = ...
    safePath.amplitudePreservationError;
report.safePath.pathGraph = safePath.pathGraph;
report.safePath.correctedSignal = squeeze(real(sum( ...
    safePath.correctedSelectedComponent,3)));
report.safePath.pathPhasorDecimated = ...
    safePath.selectedPathPhasor(1:20:end,:,:);
report.safePath.decisionRule = safePath.decisionRule;
report.truth = truth;
if ~isempty(oraclePath)
    report.oraclePathDecimated = oraclePath(1:20:end,:,:);
else
    report.oraclePathDecimated = [];
end
report.metrics = metrics;
report.encoderPolicy = ['Only vibration columns 1:3 enter estimation. ', ...
    'Encoder A-phase is decoded after blind parameters freeze and is used ', ...
    'only for scoring.'];
save([prefix,'.mat'],'report','-v7.3');

fprintf('Three-channel phase run complete: %s\n',sourceInfo.baseName);
fprintf('  sensors: 0, 120, 240 deg; blind mesh %.4f Hz\n', ...
    blindMeshFrequencyHz);
fprintf('  vibration-selected harmonic orders: ');
fprintf('%d ',harmonicOrder);
fprintf('\n');
fprintf('  selected non-3k penalty: %.3g\n',selectedOffGridPenalty);
disp('  selected sensors for common phase:');
disp(selectedCommonSensorMask);
disp('  path sensor-harmonic reliability weights:');
disp(adaptiveOptions.pathSensorHarmonicWeight);
disp('  active sensor mask in final graph:');
disp(adaptiveThree.pathGraph.activeSensorMask);
fprintf('  safe path selection: %s, sensors %s\n', ...
    safePath.selectedCandidate,strjoin(string( ...
    safePath.selectedSensorIndex),'+'));
fprintf('  held-out safeguard; full coherence %.4f -> %.4f\n', ...
    safePath.fullRawCoherence,safePath.fullCorrectedCoherence);
disp(safePath.cvTable);
disp(metrics);
end

function [selected,tableOut] = selectReliableHarmonics( ...
    spectrum,f0,fs,nSelect)
candidate = 1:min(13,floor(0.45*fs/f0));
nSensor = size(spectrum.amplitude,2);
prominence = zeros(numel(candidate),nSensor);
for ih = 1:numel(candidate)
    center = candidate(ih)*f0;
    narrow = abs(spectrum.frequency-center)<=8;
    background = abs(spectrum.frequency-center)<=40 & ...
        abs(spectrum.frequency-center)>12;
    for is = 1:nSensor
        prominence(ih,is) = max(spectrum.amplitude(narrow,is))/ ...
            max(median(spectrum.amplitude(background,is)),eps);
    end
end
jointScore = exp(mean(log(max(prominence,eps)),2));
rankScore = jointScore;
rankScore(candidate==1) = inf;
[~,rank] = sort(rankScore,'descend');
selected = sort(candidate(rank(1:min(nSelect,numel(rank)))));
isSelected = ismember(candidate(:),selected(:));
variableName = arrayfun(@(x)sprintf('prominence_c%d',x),1:nSensor, ...
    'UniformOutput',false);
tableOut = array2table([candidate(:),prominence,jointScore,isSelected], ...
    'VariableNames',[{'harmonic_order'},variableName, ...
    {'geometric_mean_prominence','selected'}]);
end

function compact = compactPhaseResult(input)
compact.meshPhase = input.meshPhase;
compact.carrierPhase = input.carrierPhase;
compact.instantaneousMeshFrequencyHz = input.instantaneousMeshFrequencyHz;
compact.pathGraph = input.pathGraph;
compact.spatialOrder = input.spatialOrder;
compact.pathPenalty = input.pathPenalty;
compact.options = input.options;
compact.identifiability = input.identifiability;
end

function reliability = sensorReliabilityFromGraph(graph,nSensor,nHarmonic)
reliability = zeros(nSensor,nHarmonic);
for ih = 1:nHarmonic
    for is = 1:nSensor
        use = graph.pairIndex(:,1)==is | graph.pairIndex(:,2)==is;
        reliability(is,ih) = max(graph.edgeFitCoherence(use,ih));
    end
end
reliability(reliability<0.50) = 0;
reliability = reliability.^2;
for ih = 1:nHarmonic
    if any(reliability(:,ih)>0)
        reliability(:,ih) = reliability(:,ih)/max(reliability(:,ih));
    end
end
end

function [selectedMask,tableOut] = selectCommonSensorSubset( ...
    component,fs,harmonicOrder,f0,smoothSec)
nSensor = size(component,2);
if nSensor==2
    subset = true(1,2);
else
    subset = logical([1 1 0;1 0 1;0 1 1;1 1 1]);
end
nSubset = size(subset,1);
frequencyRippleRmsHz = zeros(nSubset,1);
medianFrequencyHz = zeros(nSubset,1);
selectionScore = zeros(nSubset,1);
for ic = 1:nSubset
    phase = blindCommonPhase(component(:,subset(ic,:),:),fs, ...
        harmonicOrder,f0,smoothSec);
    frequency = instantaneousFrequency(phase,fs,0.10);
    edge = round(0.50*fs);
    use = (1+edge):(numel(frequency)-edge);
    frequency = frequency(use);
    trend = movingAverage(frequency,round(0.50*fs));
    frequencyRippleRmsHz(ic) = sqrt(mean((frequency-trend).^2));
    medianFrequencyHz(ic) = median(frequency);
    selectionScore(ic) = sqrt(frequencyRippleRmsHz(ic)^2+ ...
        (medianFrequencyHz(ic)-f0)^2)/sqrt(sum(subset(ic,:)));
end
[~,best] = min(selectionScore);
selectedMask = subset(best,:);
sensorSubset = strings(nSubset,1);
for ic = 1:nSubset
    sensorSubset(ic) = strjoin(string(find(subset(ic,:))),'+');
end
selected = false(nSubset,1);
selected(best) = true;
tableOut = table(sensorSubset,frequencyRippleRmsHz,medianFrequencyHz, ...
    selectionScore,selected);
end

function [signal,encoderA,fs,index,info] = ...
    loadThreeChannelSegment(file,durationSec)
[~,baseName,extension] = fileparts(file);
info.baseName = [baseName,extension];
info.nominalRpm = parseNominalRpm(info.baseName);
loaded = load(file);
if isfield(loaded,'chanvals')
    fs = 51200;
    nSample = size(loaded.chanvals,1);
    startSec = 30;
    first = min(round(startSec*fs)+1,nSample-round(durationSec*fs));
    last = min(nSample,first+round(durationSec*fs)-1);
    index = (first:last).';
    signal = loaded.chanvals(index,1:3);
    encoderA = loaded.chanvals(:,4);
    info.layout = 'chanvals: vibration 1:3, encoder A/B 4:5';
elseif isfield(loaded,'Data')
    if isfield(loaded,'SampleRate')
        fs = double(loaded.SampleRate);
    else
        fs = 100000;
    end
    nSample = size(loaded.Data,1);
    startSec = 5;
    first = min(round(startSec*fs)+1,nSample-round(durationSec*fs));
    last = min(nSample,first+round(durationSec*fs)-1);
    index = (first:last).';
    signal = loaded.Data(index,1:3);
    encoderA = [];
    info.layout = 'Data: vibration 1:3; remaining channels not used';
else
    error('run_real_three_channel_phase:Variable', ...
        'Expected chanvals or Data in %s.',file);
end
info.fs = fs;
info.segmentStartSec = (index(1)-1)/fs;
info.segmentDurationSec = numel(index)/fs;
clear loaded;
end

function rpm = parseNominalRpm(name)
token = regexp(name,'(?:S|_)(300|600|900)(?:RPM|L|_)', ...
    'tokens','once');
if isempty(token)
    rpm = 600;
else
    rpm = str2double(token{1});
end
end

function [frequencyHz,out] = blindMeshPeak(signal,fs,searchBand)
n = size(signal,1);
window = hann(n);
spectrum = 2*abs(fft(signal.*window,[],1))/sum(window);
frequency = (0:n-1)'*fs/n;
amplitude = sqrt(mean(spectrum.^2,2));
use = frequency>=searchBand(1) & frequency<=searchBand(2);
candidateFrequency = frequency(use);
candidateAmplitude = amplitude(use);
[~,peak] = max(candidateAmplitude);
frequencyHz = candidateFrequency(peak);
keep = 1:floor(n/2);
out.frequency = frequency(keep);
out.amplitude = spectrum(keep,:);
out.combinedAmplitude = amplitude(keep);
end

function phase = blindCommonPhase(component,fs,harmonicOrder,f0,smoothSec)
if ismatrix(component)
    component = reshape(component,size(component,1),size(component,2),1);
end
[nSample,nSensor,nHarmonic] = size(component);
nominal = 2*pi*f0*(0:nSample-1)'/fs;
numerator = zeros(nSample-1,1);
denominator = zeros(nSample-1,1);
for is = 1:nSensor
    for ih = 1:nHarmonic
        amplitude = abs(component(:,is,ih));
        normalized = amplitude/max(median(amplitude),eps);
        weight = normalized.^2./(1+normalized.^2);
        pairWeight = min(weight(2:end),weight(1:end-1));
        increment = angle(component(2:end,is,ih).* ...
            conj(component(1:end-1,is,ih)).* ...
            exp(-1i*harmonicOrder(ih)*diff(nominal)))/harmonicOrder(ih);
        numerator = numerator+pairWeight.*increment;
        denominator = denominator+pairWeight;
    end
end
increment = numerator./max(denominator,eps);
increment = movingAverage(increment,round(0.0015*fs));
raw = nominal+[0;cumsum(increment)];
phase = nominal+movingAverage(raw-nominal,round(smoothSec*fs));
end

function [selected,tableOut] = selectOffGridPenaltyMulti(component, ...
    carrierPhase,maxOrder,spatialMultiple,ridge,candidate,edgeReliability)
stride = 20;
index = (1:stride:numel(carrierPhase)).';
component = component(index,:,:);
carrierPhase = carrierPhase(index);
order = -maxOrder:maxOrder;
basis = exp(1i*carrierPhase*order);
cycle = floor((carrierPhase-carrierPhase(1))/(2*pi));
train = mod(cycle,2)==0;
test = ~train;
[~,nSensor,nHarmonic] = size(component);
pair = nchoosek(1:nSensor,2);
validationRmseDeg = zeros(numel(candidate),1);
for ic = 1:numel(candidate)
    penalty = ones(size(order));
    penalty(mod(abs(order),spatialMultiple)~=0) = candidate(ic);
    sumSquare = 0;
    sumWeight = 0;
    for ih = 1:nHarmonic
        for ip = 1:size(pair,1)
            sensorI = pair(ip,1);
            sensorJ = pair(ip,2);
            cross = component(:,sensorJ,ih).*conj(component(:,sensorI,ih));
            unitCross = cross./max(abs(cross),eps);
            amplitude = sqrt(abs(component(:,sensorI,ih)).* ...
                abs(component(:,sensorJ,ih)));
            normalized = amplitude/max(median(amplitude),eps);
            weight = normalized.^2./(1+normalized.^2)* ...
                max(edgeReliability(ip,ih)^2,0.01);
            coefficient = complexFit(basis(train,:),unitCross(train), ...
                weight(train),ridge,penalty);
            prediction = basis(test,:)*coefficient;
            prediction = prediction./max(abs(prediction),eps);
            error = angle(unitCross(test).*conj(prediction));
            sumSquare = sumSquare+sum(weight(test).*error.^2);
            sumWeight = sumWeight+sum(weight(test));
        end
    end
    validationRmseDeg(ic) = rad2deg(sqrt(sumSquare/sumWeight));
end
[~,best] = min(validationRmseDeg);
selected = candidate(best);
tableOut = table(candidate(:),validationRmseDeg, ...
    candidate(:)==selected,'VariableNames',{'off_grid_penalty', ...
    'validation_phase_rmse_deg','selected'});
end

function path = fitKnownCarrierPath(component,carrierPhase,order,ridge)
[nSample,nSensor,nHarmonic] = size(component);
pair = nchoosek(1:nSensor,2);
nPair = size(pair,1);
incidence = zeros(nPair,nSensor);
for ip = 1:nPair
    incidence(ip,pair(ip,1)) = -1;
    incidence(ip,pair(ip,2)) = 1;
end
basis = exp(1i*carrierPhase*order);
path = complex(ones(size(component)));
for ih = 1:nHarmonic
    edgePhase = zeros(nSample,nPair);
    edgeWeight = zeros(nPair,1);
    for ip = 1:nPair
        sensorI = pair(ip,1);
        sensorJ = pair(ip,2);
        cross = component(:,sensorJ,ih).*conj(component(:,sensorI,ih));
        unitCross = cross./max(abs(cross),eps);
        amplitude = sqrt(abs(component(:,sensorI,ih)).* ...
            abs(component(:,sensorJ,ih)));
        normalized = amplitude/max(median(amplitude),eps);
        weight = normalized.^2./(1+normalized.^2);
        coefficient = complexFit(basis,unitCross,weight,ridge, ...
            ones(size(order)));
        fitted = basis*coefficient;
        fitted = fitted./max(abs(fitted),eps);
        edgePhase(:,ip) = unwrap(angle(fitted));
        error = angle(unitCross.*conj(fitted));
        edgeWeight(ip) = 1/max(sum(weight.*error.^2)/sum(weight), ...
            deg2rad(1)^2);
    end
    edgePhase = alignPairBranches(edgePhase,pair,nSensor);
    edgeWeight = min(max(edgeWeight/median(edgeWeight),0.1),10);
    laplacian = incidence'*diag(edgeWeight)*incidence;
    graphMap = (laplacian+ones(nSensor)/nSensor)\ ...
        (incidence'*diag(edgeWeight));
    sensorPhase = edgePhase*graphMap.';
    for is = 1:nSensor
        path(:,is,ih) = exp(1i*sensorPhase(:,is));
    end
end
end

function edgePhase = alignPairBranches(edgePhase,pair,nSensor)
reference = zeros(size(edgePhase,1),nSensor);
for sensorJ = 2:nSensor
    edge = find(pair(:,1)==1 & pair(:,2)==sensorJ,1);
    reference(:,sensorJ) = edgePhase(:,edge);
end
for ip = 1:size(pair,1)
    target = reference(:,pair(ip,2))-reference(:,pair(ip,1));
    edgePhase(:,ip) = edgePhase(:,ip)+2*pi*round( ...
        median((target-edgePhase(:,ip))/(2*pi)));
end
end

function coefficient = complexFit(basis,target,weight,ridge,penalty)
rootWeight = sqrt(max(weight,0));
weightedBasis = basis.*rootWeight;
normal = weightedBasis'*weightedBasis;
coefficient = (normal+ridge*trace(normal)/size(normal,1)* ...
    diag(penalty))\(weightedBasis'*(target.*rootWeight));
end

function [rmseDeg,p95Deg,carrierRmseRpm] = phaseScore( ...
    estimate,truth,fs,zr)
error = estimate-truth.meshPhase;
error = error-median(error);
errorDeg = rad2deg(error);
rmseDeg = sqrt(mean(errorDeg.^2));
p95Deg = percentile(abs(errorDeg),95);
estimatedFrequency = instantaneousFrequency(estimate,fs,0.10);
edge = round(0.20*fs);
use = (1+edge):(numel(estimate)-edge);
frequencyRmse = sqrt(mean((estimatedFrequency(use)- ...
    truth.meshFrequencyHz(use)).^2));
carrierRmseRpm = frequencyRmse/zr*60;
end

function value = pairResidual(component)
[~,nSensor,nHarmonic] = size(component);
pair = nchoosek(1:nSensor,2);
residual = zeros(size(pair,1),nHarmonic);
for ip = 1:size(pair,1)
    for ih = 1:nHarmonic
        difference = angle(component(:,pair(ip,2),ih).* ...
            conj(component(:,pair(ip,1),ih)));
        center = angle(mean(exp(1i*difference)));
        residual(ip,ih) = sqrt(mean(angle(exp(1i* ...
            (difference-center))).^2));
    end
end
value = rad2deg(mean(residual,'all'));
end

function value = pairCoherence(component)
[~,nSensor,nHarmonic] = size(component);
pair = nchoosek(1:nSensor,2);
coherence = zeros(size(pair,1),nHarmonic);
for ip = 1:size(pair,1)
    for ih = 1:nHarmonic
        sensorI = pair(ip,1);
        sensorJ = pair(ip,2);
        cross = component(:,sensorJ,ih).*conj(component(:,sensorI,ih));
        weight = sqrt(abs(component(:,sensorJ,ih)).* ...
            abs(component(:,sensorI,ih)));
        coherence(ip,ih) = abs(sum(weight.*cross./ ...
            max(abs(cross),eps))/sum(weight));
    end
end
value = mean(coherence,'all');
end

function rmseDeg = pathScore(estimatePath,truePath,component,graph)
[~,nSensor,nHarmonic] = size(component);
pair = nchoosek(1:nSensor,2);
sumSquare = 0;
sumWeight = 0;
for ip = 1:size(pair,1)
    for ih = 1:nHarmonic
        if ~isempty(graph) && ~graph.activeEdgeMask(ip,ih)
            continue;
        end
        sensorI = pair(ip,1);
        sensorJ = pair(ip,2);
        estimate = estimatePath(:,sensorJ,ih).* ...
            conj(estimatePath(:,sensorI,ih));
        truth = truePath(:,sensorJ,ih).*conj(truePath(:,sensorI,ih));
        error = angle(estimate.*conj(truth));
        weight = sqrt(abs(component(:,sensorI,ih)).* ...
            abs(component(:,sensorJ,ih)));
        sumSquare = sumSquare+sum(weight.*error.^2);
        sumWeight = sumWeight+sum(weight);
    end
end
if sumWeight==0
    rmseDeg = nan;
else
    rmseDeg = rad2deg(sqrt(sumSquare/sumWeight));
end
end

function value = graphPairCoherence(component,graph)
[~,~,nHarmonic] = size(component);
sumValue = 0;
count = 0;
for ip = 1:size(graph.pairIndex,1)
    for ih = 1:nHarmonic
        if ~graph.activeEdgeMask(ip,ih)
            continue;
        end
        sensorI = graph.pairIndex(ip,1);
        sensorJ = graph.pairIndex(ip,2);
        cross = component(:,sensorJ,ih).*conj(component(:,sensorI,ih));
        weight = sqrt(abs(component(:,sensorJ,ih)).* ...
            abs(component(:,sensorI,ih)));
        sumValue = sumValue+abs(sum(weight.*cross./ ...
            max(abs(cross),eps))/sum(weight));
        count = count+1;
    end
end
if count==0
    value = nan;
else
    value = sumValue/count;
end
end

function edge = risingEdges(signal,threshold,minSeparation)
edge = find(diff(signal>threshold)==1)+1;
edge = edge([true;diff(edge)>minSeparation]);
end

function frequency = instantaneousFrequency(phase,fs,smoothSec)
frequency = [diff(phase);phase(end)-phase(end-1)]*fs/(2*pi);
frequency = movingAverage(frequency,round(smoothSec*fs));
end

function y = movingAverage(x,window)
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

function [centerDeg,phaseDeg] = binnedPhaseMap(carrierPhase,phasor,nBin)
wrapped = mod(carrierPhase,2*pi);
edge = linspace(0,2*pi,nBin+1);
centerDeg = rad2deg((edge(1:end-1)+edge(2:end))/2).';
phaseDeg = nan(nBin,1);
for ib = 1:nBin
    use = wrapped>=edge(ib) & wrapped<edge(ib+1);
    if any(use)
        phaseDeg(ib) = rad2deg(angle(mean(phasor(use))));
    end
end
end
