function report = run_real_tacholess_joint_phase(tdmsFile)
%RUN_REAL_TACHOLESS_JOINT_PHASE Blind phase separation on measured data.
% Vibration columns 3/4 are the only estimator inputs. Tacho columns 1/2
% are isolated until the final scoring block.

setup_paths;
rootDir = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || strlength(string(tdmsFile))==0
    tdmsFile = fullfile(rootDir,'data_local','sun','dual_channel_sun50.tdms');
end
fs = 51200;
zs = 21;
zr = 84;
sunEncoderPpr = 120;
segmentStartSec = 3.0;
segmentDurationSec = 12.0;

data = tdmsread(tdmsFile);
if numel(data)~=1 || width(data{1})<4
    error('run_real_tacholess_joint_phase:ChannelLayout', ...
        'Expected one TDMS group with two tacho and two vibration channels.');
end
tableData = data{1};
carrierIndex = table2array(tableData(:,1));
sunTacho = table2array(tableData(:,2));
vibration = [table2array(tableData(:,3)),table2array(tableData(:,4))];
vibration = vibration-mean(vibration,1);

first = round(segmentStartSec*fs)+1;
last = min(size(vibration,1),first+round(segmentDurationSec*fs)-1);
segmentIndex = (first:last).';
rawSignal = vibration(segmentIndex,:);
tSegment = (0:numel(segmentIndex)-1)'/fs;

% Vibration-only nominal mesh-frequency search. No pulse channel is touched.
[blindMeshFrequencyHz,blindSpectrum] = blindMeshPeak(rawSignal,fs,[150,185]);
harmonicOrder = [1,4,11,13];
centerFrequencyHz = harmonicOrder*blindMeshFrequencyHz;
halfBandwidthHz = [36,52,82,92];
prominence = harmonicProminence(blindSpectrum,harmonicOrder, ...
    blindMeshFrequencyHz);

extracted = pg_extract_analytic_harmonics(rawSignal,fs, ...
    centerFrequencyHz,halfBandwidthHz);
edgeSec = 0.50;
use = tSegment>=edgeSec & tSegment<=tSegment(end)-edgeSec;
t = tSegment(use)-tSegment(find(use,1));
component = extracted.component(use,:,:);
signalUse = rawSignal(use,:);

options.nIteration = 3;
options.spatialMultiple = 3;
options.pathOrder = 4;
options.ridge = 2e-3;
options.relaxation = 0.35;
options.commonAnchorWeight = 0.80;
options.phaseSmoothSec = 0.012;
options.frequencySmoothSec = 0.10;
options.removeCommonSpatialOrder = false;
multiNoPathPhase = commonPhaseIncrement(component,fs,harmonicOrder, ...
    blindMeshFrequencyHz,0.012);
singleH1Phase = commonPhaseIncrement(component(:,1,1),fs,1, ...
    blindMeshFrequencyHz,0.012);

frequencyCorrectionCandidateHz = -0.12:0.002:0.12;
[selectedFrequencyCorrectionHz,frequencyCorrectionCvTable, ...
    frequencyCorrectionAccepted] = ...
    selectFrequencyCorrection(component,multiNoPathPhase,t,zr,3,4, ...
    2e-3,frequencyCorrectionCandidateHz);
options.meanFrequencyCorrectionHz = selectedFrequencyCorrectionHz;
joint = pg_joint_tacholess_phase(component,fs,harmonicOrder, ...
    blindMeshFrequencyHz,zr,options);
writetable(frequencyCorrectionCvTable,fullfile('results', ...
    'real_path_repeatability_frequency_cv.csv'));

freeOptions = options;
freeOptions.spatialMultiple = 1;
freeOptions.pathOrder = 12;
unconstrained = pg_joint_tacholess_phase(component,fs,harmonicOrder, ...
    blindMeshFrequencyHz,zr,freeOptions);

% Soft three-planet prior: select non-3k shrinkage by odd/even carrier-cycle
% cross-validation using vibration only. This permits repeatable casing and
% mounting asymmetry without discarding the geometric 3k core.
penaltyCandidate = [1,3,10,30,100,300];
[selectedOffGridPenalty,penaltyCvTable] = selectOffGridPenalty( ...
    component,(multiNoPathPhase+2*pi*selectedFrequencyCorrectionHz*t)/zr, ...
    12,3,2e-3,penaltyCandidate);
hybridOptions = options;
hybridOptions.spatialOrderOverride = -12:12;
hybridOptions.offGridPenalty = selectedOffGridPenalty;
hybrid = pg_joint_tacholess_phase(component,fs,harmonicOrder, ...
    blindMeshFrequencyHz,zr,hybridOptions);
writetable(penaltyCvTable,fullfile('results', ...
    'real_hybrid_path_penalty_cv.csv'));

%% Encoder truth: this block is intentionally after every blind estimate.
sunEdges = risingEdges(sunTacho,1.0,5);
sunPulsePhase = (0:numel(sunEdges)-1)'*(2*pi/sunEncoderPpr);
sampleNumber = (1:size(vibration,1)).';
sunPhase = interp1(sunEdges,sunPulsePhase,sampleNumber,'linear','extrap');
meshPerSun = zs*zr/(zs+zr);
encoderMeshPhase = meshPerSun*sunPhase;
encoderCarrierPhase = encoderMeshPhase/zr;

truth.meshPhase = encoderMeshPhase(segmentIndex(use));
truth.carrierPhase = encoderCarrierPhase(segmentIndex(use));
truth.meshFrequencyHz = instantaneousFrequency( ...
    truth.meshPhase,fs,0.10);
truth.sunEdges = sunEdges;
truth.definition = ['Column 2, 120 PPR sun encoder; fixed-ring kinematics ', ...
    'phi_m = Zs*Zr/(Zs+Zr)*phi_s. Not supplied to the estimator.'];

carrierEdges = risingEdges(carrierIndex,2.0,100);
carrierRevolution = diff(encoderCarrierPhase(carrierEdges))/(2*pi);
encoderConsistency.revolutionMedian = median(carrierRevolution);
encoderConsistency.revolutionRmse = ...
    sqrt(mean((carrierRevolution-1).^2));
encoderConsistency.carrierEdgeCount = numel(carrierEdges);

% Encoder-angle path maps are oracles used only to score the frozen blind
% maps. The full -12:12 oracle is the common comparison target; the hard-3k
% oracle is retained only to diagnose the geometric subspace itself.
oraclePathHard = fitRelativePath(component,truth.carrierPhase,3, ...
    options.pathOrder,2e-3);
oraclePathFull = fitRelativePath(component,truth.carrierPhase,1, ...
    freeOptions.pathOrder,2e-3);

method = {'single_channel_H1';'multiharmonic_no_path'; ...
    'unconstrained_joint';'hard_three_planet_joint'; ...
    'adaptive_three_planet_joint'};
phaseSet = {singleH1Phase,multiNoPathPhase, ...
    unconstrained.meshPhase,joint.meshPhase,hybrid.meshPhase};
componentSet = {component,component,unconstrained.correctedComponent, ...
    joint.correctedComponent,hybrid.correctedComponent};
pathSet = {[],[],unconstrained.pathPhasor,joint.pathPhasor, ...
    hybrid.pathPhasor};
parameterCount = [0;0;numel(harmonicOrder)*(2*freeOptions.pathOrder+1); ...
    numel(harmonicOrder)*(2*options.pathOrder+1); ...
    numel(harmonicOrder)*numel(hybridOptions.spatialOrderOverride)];

nMethod = numel(method);
meshPhaseRmseDeg = zeros(nMethod,1);
meshPhaseP95Deg = zeros(nMethod,1);
meshFrequencyRmseHz = zeros(nMethod,1);
carrierSpeedRmseRpm = zeros(nMethod,1);
dualPhaseResidualDeg = zeros(nMethod,1);
dualPhaseCoherence = zeros(nMethod,1);
pathVsEncoderRmseDeg = nan(nMethod,1);
amplitudePreservationError = zeros(nMethod,1);
for im = 1:nMethod
    [meshPhaseRmseDeg(im),meshPhaseP95Deg(im), ...
        meshFrequencyRmseHz(im),carrierSpeedRmseRpm(im)] = ...
        phaseScore(phaseSet{im},truth,fs,zr);
    dualPhaseResidualDeg(im) = dualResidual(componentSet{im});
    dualPhaseCoherence(im) = dualCoherence(componentSet{im});
    amplitudePreservationError(im) = max(abs(abs(componentSet{im}(:))- ...
        abs(component(:))))/max(abs(component(:)));
    if ~isempty(pathSet{im})
        pathVsEncoderRmseDeg(im) = pathDifferenceScore( ...
            pathSet{im},oraclePathFull,component);
    end
end

metrics = table(method,parameterCount,meshPhaseRmseDeg,meshPhaseP95Deg, ...
    meshFrequencyRmseHz,carrierSpeedRmseRpm,pathVsEncoderRmseDeg, ...
    dualPhaseResidualDeg,dualPhaseCoherence,amplitudePreservationError, ...
    'VariableNames',{'method','complex_path_parameter_count', ...
    'mesh_phase_rmse_deg','mesh_phase_p95_deg', ...
    'mesh_frequency_rmse_hz','carrier_speed_rmse_rpm', ...
    'path_phase_vs_encoder_rmse_deg','dual_phase_residual_deg', ...
    'dual_phase_coherence','amplitude_preservation_error'});
writetable(metrics,fullfile('results','real_tacholess_joint_phase_metrics.csv'));

% Test whether measured differential paths actually occupy 3k orders.
[spatialOrderTable,threePlanetEnergyFraction] = spatialOrderAudit( ...
    component,truth.carrierPhase,12,2e-3,harmonicOrder);
writetable(spatialOrderTable,fullfile('results', ...
    'real_three_planet_spatial_orders.csv'));

rawCarrierMatrix = pg_phase_synchronous_matrix(signalUse, ...
    hybrid.carrierPhase,2048);
correctedCarrierMatrix = pg_phase_synchronous_matrix( ...
    hybrid.correctedSignal,hybrid.carrierPhase,2048);

fig = figure('Visible','off','Color','w','Position',[50 40 1220 830]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
plot(t,movingAverage(truth.meshFrequencyHz,round(0.50*fs)), ...
    'k','LineWidth',1.2); hold on;
plot(t,movingAverage(hybrid.instantaneousMeshFrequencyHz,round(0.50*fs)), ...
    'Color',[0.85 0.20 0.10], ...
    'LineWidth',0.9);
xlim([t(1),t(end)]); grid on; box off;
xlabel('Time (s)'); ylabel('Mesh frequency (Hz)');
title('Measured data: blind speed recovery');
legend('Hidden encoder truth','Adaptive three-planet joint','Location','best');

nexttile;
displayMethod = [2,3,4,5];
color = lines(nMethod);
for im = displayMethod
    plot(t,alignedPhaseErrorDeg(phaseSet{im},truth.meshPhase), ...
        'Color',color(im,:),'LineWidth',0.8); hold on;
end
xlim([t(1),t(end)]); grid on; box off;
xlabel('Time (s)'); ylabel('Mesh-phase error (deg)');
title('Encoder-hidden phase error');
legend(method(displayMethod),'Interpreter','none','Location','best');

nexttile;
ihPlot = find(harmonicOrder==11,1);
oracleDifference = angle(oraclePathFull(:,2,ihPlot).* ...
    conj(oraclePathFull(:,1,ihPlot)));
blindDifference = angle(hybrid.pathPhasor(:,2,ihPlot).* ...
    conj(hybrid.pathPhasor(:,1,ihPlot)));
[carrierAngleDeg,oracleMapDeg] = binnedPhaseMap( ...
    truth.carrierPhase,exp(1i*oracleDifference),180);
[~,blindMapDeg] = binnedPhaseMap( ...
    truth.carrierPhase,exp(1i*blindDifference),180);
plot(carrierAngleDeg,oracleMapDeg,'k','LineWidth',1.1); hold on;
plot(carrierAngleDeg,blindMapDeg,'Color',[0.85 0.20 0.10], ...
    'LineWidth',0.8);
xlim([0 360]); ylim([-180 180]); grid on; box off;
xlabel('Encoder carrier angle for scoring (deg)');
ylabel('90-0 deg path phase (deg)');
title('Moving-path phase: encoder oracle vs blind estimate');
legend('Encoder-angle oracle','Blind estimate','Location','best');

nexttile;
figureMethod = {'H1','Multi-H','Free','Hard 3k','Adaptive 3k'};
methodCategory = categorical(figureMethod,figureMethod,'Ordinal',true);
yyaxis left;
bar(methodCategory,dualPhaseResidualDeg,'FaceColor',[0.25 0.55 0.80]);
ylabel('Dual-channel residual (deg)');
yyaxis right;
plot(methodCategory,dualPhaseCoherence,'ro-', ...
    'LineWidth',1.2,'MarkerFaceColor','r');
ylabel('Dual-channel phase coherence'); ylim([0 1.05]);
grid on; box off; title('Measured differential-path alignment');
exportgraphics(fig,fullfile('results', ...
    'real_tacholess_joint_phase_validation.png'),'Resolution',190);
close(fig);

fig = figure('Visible','off','Color','w','Position',[80 70 1180 720]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
order = unique(spatialOrderTable.spatial_order).';
for ih = 1:numel(harmonicOrder)
    nexttile;
    useOrder = spatialOrderTable.harmonic_order==harmonicOrder(ih);
    coefficient = spatialOrderTable.coefficient_magnitude(useOrder);
    bar(order,coefficient/max(coefficient),'FaceColor',[0.20 0.50 0.75]);
    hold on;
    threeOrder = order(mod(abs(order),3)==0);
    stem(threeOrder,ones(size(threeOrder)),'r:','Marker','none');
    xlim([order(1)-0.5,order(end)+0.5]); ylim([0 1.08]);
    grid on; box off; xlabel('Carrier spatial order');
    ylabel('Normalized coefficient');
    title(sprintf('h=%d, 3k energy fraction %.3f',harmonicOrder(ih), ...
        threePlanetEnergyFraction(ih)));
end
exportgraphics(fig,fullfile('results', ...
    'real_three_planet_spatial_orders.png'),'Resolution',190);
close(fig);

report.sourceFile = tdmsFile;
report.fs = fs;
report.segmentStartSec = segmentStartSec;
report.segmentDurationSec = segmentDurationSec;
report.t = t;
report.rawSignal = signalUse;
report.blindMeshFrequencyHz = blindMeshFrequencyHz;
report.harmonicOrder = harmonicOrder;
report.harmonicProminence = prominence;
report.truth = truth;
report.encoderConsistency = encoderConsistency;
report.joint = joint;
report.unconstrained = unconstrained;
report.hybrid = hybrid;
report.selectedFrequencyCorrectionHz = selectedFrequencyCorrectionHz;
report.frequencyCorrectionAccepted = frequencyCorrectionAccepted;
report.frequencyCorrectionCvTable = frequencyCorrectionCvTable;
report.selectedOffGridPenalty = selectedOffGridPenalty;
report.penaltyCvTable = penaltyCvTable;
report.oraclePathFull = oraclePathFull;
report.oraclePathHard = oraclePathHard;
report.metrics = metrics;
report.spatialOrderTable = spatialOrderTable;
report.threePlanetEnergyFraction = threePlanetEnergyFraction;
report.rawCarrierMatrix = rawCarrierMatrix;
report.correctedCarrierMatrix = correctedCarrierMatrix;
report.encoderPolicy = ['TDMS columns 1/2 are read only after all blind ', ...
    'estimates; they are used for scoring and spatial-order audit only.'];
save(fullfile('results','real_tacholess_joint_phase.mat'),'report','-v7.3');

encoderMeanMeshHz = mean(truth.meshFrequencyHz);
fprintf('Real tacholess phase run complete.\n');
fprintf('  vibration-only mesh estimate %.4f Hz; encoder mean %.4f Hz\n', ...
    blindMeshFrequencyHz,encoderMeanMeshHz);
fprintf('  carrier-index consistency median %.6f rev, RMSE %.6g rev\n', ...
    encoderConsistency.revolutionMedian, ...
    encoderConsistency.revolutionRmse);
fprintf('  three-planet 3k energy fractions: ');
fprintf('%.3f ',threePlanetEnergyFraction);
fprintf('\n');
fprintf('  vibration-only selected off-grid penalty: %.3g\n', ...
    selectedOffGridPenalty);
fprintf('  path-repeatability frequency correction: %+.6f Hz\n', ...
    selectedFrequencyCorrectionHz);
fprintf('  frequency correction accepted: %d\n',frequencyCorrectionAccepted);
disp(metrics);
end

function [meshFrequencyHz,out] = blindMeshPeak(signal,fs,searchBand)
n = size(signal,1);
window = hann(n);
spectrum = 2*abs(fft(signal.*window,[],1))/sum(window);
frequency = (0:n-1)'*fs/n;
amplitude = sqrt(mean(spectrum.^2,2));
use = frequency>=searchBand(1) & frequency<=searchBand(2);
candidateFrequency = frequency(use);
candidateAmplitude = amplitude(use);
[~,index] = max(candidateAmplitude);
meshFrequencyHz = candidateFrequency(index);
out.frequency = frequency(1:floor(n/2));
out.amplitude = spectrum(1:floor(n/2),:);
out.combinedAmplitude = amplitude(1:floor(n/2));
out.searchBandHz = searchBand;
end

function tableOut = harmonicProminence(spectrum,harmonicOrder,f0)
nHarmonic = numel(harmonicOrder);
peakFrequencyHz = zeros(nHarmonic,2);
prominence = zeros(nHarmonic,2);
for ih = 1:nHarmonic
    center = harmonicOrder(ih)*f0;
    narrow = abs(spectrum.frequency-center)<=12;
    background = abs(spectrum.frequency-center)<=60 & ~narrow;
    for is = 1:2
        [peak,index] = max(spectrum.amplitude(narrow,is));
        localFrequency = spectrum.frequency(narrow);
        peakFrequencyHz(ih,is) = localFrequency(index);
        prominence(ih,is) = peak/median(spectrum.amplitude(background,is));
    end
end
tableOut = table(harmonicOrder(:),peakFrequencyHz(:,1), ...
    peakFrequencyHz(:,2),prominence(:,1),prominence(:,2), ...
    'VariableNames',{'harmonic_order','peak_0_hz','peak_90_hz', ...
    'prominence_0','prominence_90'});
end

function edge = risingEdges(signal,threshold,minSeparation)
edge = find(diff(signal>threshold)==1)+1;
edge = edge([true;diff(edge)>minSeparation]);
end

function phase = commonPhaseIncrement(component,fs,harmonicOrder,f0,smoothSec)
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

function path = fitRelativePath(component,carrierPhase,spatialMultiple, ...
    pathOrder,ridge)
[nSample,nSensor,nHarmonic] = size(component);
order = spatialMultiple*(-pathOrder:pathOrder);
basis = exp(1i*carrierPhase*order);
path = complex(ones(size(component)));
for ih = 1:nHarmonic
    relativePhase = zeros(nSample,nSensor);
    for is = 2:nSensor
        cross = component(:,is,ih).*conj(component(:,1,ih));
        unitCross = cross./max(abs(cross),eps);
        amplitude = sqrt(abs(component(:,is,ih)).*abs(component(:,1,ih)));
        normalized = amplitude/max(median(amplitude),eps);
        weight = normalized.^2./(1+normalized.^2);
        coefficient = complexFit(basis,unitCross,weight,ridge);
        fitted = basis*coefficient;
        relativePhase(:,is) = unwrap(angle(fitted));
    end
    meanRelative = mean(relativePhase,2);
    for is = 1:nSensor
        path(:,is,ih) = exp(1i*(relativePhase(:,is)-meanRelative));
    end
end
end

function [tableOut,energyFraction] = spatialOrderAudit(component, ...
    carrierPhase,maxOrder,ridge,harmonicOrder)
order = -maxOrder:maxOrder;
basis = exp(1i*carrierPhase*order);
nHarmonic = numel(harmonicOrder);
coefficientMagnitude = zeros(nHarmonic,numel(order));
energyFraction = zeros(nHarmonic,1);
for ih = 1:nHarmonic
    cross = component(:,2,ih).*conj(component(:,1,ih));
    unitCross = cross./max(abs(cross),eps);
    amplitude = sqrt(abs(component(:,2,ih)).*abs(component(:,1,ih)));
    normalized = amplitude/max(median(amplitude),eps);
    weight = normalized.^2./(1+normalized.^2);
    coefficient = complexFit(basis,unitCross,weight,ridge);
    coefficientMagnitude(ih,:) = abs(coefficient);
    isThree = mod(abs(order),3)==0;
    energyFraction(ih) = sum(abs(coefficient(isThree)).^2)/ ...
        sum(abs(coefficient).^2);
end
tableOut = table(repelem(harmonicOrder(:),numel(order)), ...
    repmat(order(:),nHarmonic,1),reshape(coefficientMagnitude.',[],1), ...
    'VariableNames',{'harmonic_order','spatial_order', ...
    'coefficient_magnitude'});
end

function coefficient = complexFit(basis,target,weight,ridge)
rootWeight = sqrt(max(weight,0));
weightedBasis = basis.*rootWeight;
normal = weightedBasis'*weightedBasis;
coefficient = (normal+ridge*trace(normal)/size(normal,1)* ...
    eye(size(normal)))\(weightedBasis'*(target.*rootWeight));
end

function [selected,tableOut] = selectOffGridPenalty(component, ...
    carrierPhase,maxOrder,spatialMultiple,ridge,candidate)
order = -maxOrder:maxOrder;
basis = exp(1i*carrierPhase*order);
cycle = floor((carrierPhase-carrierPhase(1))/(2*pi));
train = mod(cycle,2)==0;
test = ~train;
nHarmonic = size(component,3);
validationRmseDeg = zeros(numel(candidate),1);
for ic = 1:numel(candidate)
    penalty = ones(size(order));
    penalty(mod(abs(order),spatialMultiple)~=0) = candidate(ic);
    sumSquare = 0;
    sumWeight = 0;
    for ih = 1:nHarmonic
        cross = component(:,2,ih).*conj(component(:,1,ih));
        unitCross = cross./max(abs(cross),eps);
        amplitude = sqrt(abs(component(:,2,ih)).*abs(component(:,1,ih)));
        normalized = amplitude/max(median(amplitude),eps);
        weight = normalized.^2./(1+normalized.^2);
        coefficient = complexFitWithPenalty(basis(train,:), ...
            unitCross(train),weight(train),ridge,penalty);
        prediction = basis(test,:)*coefficient;
        prediction = prediction./max(abs(prediction),eps);
        error = angle(unitCross(test).*conj(prediction));
        sumSquare = sumSquare+sum(weight(test).*error.^2);
        sumWeight = sumWeight+sum(weight(test));
    end
    validationRmseDeg(ic) = rad2deg(sqrt(sumSquare/sumWeight));
end
[~,index] = min(validationRmseDeg);
selected = candidate(index);
tableOut = table(candidate(:),validationRmseDeg, ...
    'VariableNames',{'off_grid_penalty','validation_phase_rmse_deg'});
end

function [selected,tableOut,accepted] = selectFrequencyCorrection(component, ...
    initialMeshPhase,t,zr,spatialMultiple,pathOrder,ridge,candidateHz)
% Select the mean-frequency gauge by out-of-time path repeatability. The
% calculation is deliberately decimated because the path varies at carrier
% orders, not at the 51.2-kHz waveform sampling rate.
stride = 20;
index = (1:stride:numel(t)).';
tSmall = t(index);
componentSmall = component(index,:,:);
nHarmonic = size(componentSmall,3);
order = spatialMultiple*(-pathOrder:pathOrder);
splitTime = median(tSmall);
train = tSmall<=splitTime;
test = ~train;
validationRmseDeg = zeros(numel(candidateHz),1);
for ic = 1:numel(candidateHz)
    phase = initialMeshPhase(index)+2*pi*candidateHz(ic)*tSmall;
    basis = exp(1i*(phase/zr)*order);
    sumSquare = 0;
    sumWeight = 0;
    for ih = 1:nHarmonic
        cross = componentSmall(:,2,ih).*conj(componentSmall(:,1,ih));
        unitCross = cross./max(abs(cross),eps);
        amplitude = sqrt(abs(componentSmall(:,2,ih)).* ...
            abs(componentSmall(:,1,ih)));
        normalized = amplitude/max(median(amplitude),eps);
        weight = normalized.^2./(1+normalized.^2);
        coefficient = complexFit(basis(train,:),unitCross(train), ...
            weight(train),ridge);
        prediction = basis(test,:)*coefficient;
        prediction = prediction./max(abs(prediction),eps);
        error = angle(unitCross(test).*conj(prediction));
        sumSquare = sumSquare+sum(weight(test).*error.^2);
        sumWeight = sumWeight+sum(weight(test));
    end
    validationRmseDeg(ic) = rad2deg(sqrt(sumSquare/sumWeight));
end
[~,bestIndex] = min(validationRmseDeg);
zeroIndex = find(abs(candidateHz)==min(abs(candidateHz)),1);
improvement = validationRmseDeg(zeroIndex)-validationRmseDeg(bestIndex);
accepted = bestIndex>1 && bestIndex<numel(candidateHz) && improvement>=0.20;
if accepted
    selected = candidateHz(bestIndex);
else
    selected = 0;
end
selectedColumn = candidateHz(:)==selected;
tableOut = table(candidateHz(:),validationRmseDeg,selectedColumn, ...
    'VariableNames',{'mesh_frequency_correction_hz', ...
    'path_repeatability_rmse_deg','selected_after_safety_check'});
end

function coefficient = complexFitWithPenalty(basis,target,weight,ridge,penalty)
rootWeight = sqrt(max(weight,0));
weightedBasis = basis.*rootWeight;
normal = weightedBasis'*weightedBasis;
coefficient = (normal+ridge*trace(normal)/size(normal,1)* ...
    diag(penalty))\(weightedBasis'*(target.*rootWeight));
end

function [rmseDeg,p95Deg,frequencyRmseHz,carrierRmseRpm] = ...
    phaseScore(estimate,truth,fs,zr)
errorDeg = alignedPhaseErrorDeg(estimate,truth.meshPhase);
rmseDeg = sqrt(mean(errorDeg.^2));
p95Deg = percentile(abs(errorDeg),95);
estimatedFrequency = instantaneousFrequency(estimate,fs,0.10);
edge = round(0.20*fs);
use = (1+edge):(numel(estimate)-edge);
frequencyRmseHz = sqrt(mean((estimatedFrequency(use)- ...
    truth.meshFrequencyHz(use)).^2));
carrierRmseRpm = frequencyRmseHz/zr*60;
end

function errorDeg = alignedPhaseErrorDeg(estimate,truth)
difference = estimate-truth;
difference = difference-median(difference);
errorDeg = rad2deg(difference);
end

function frequency = instantaneousFrequency(phase,fs,smoothSec)
frequency = [diff(phase);phase(end)-phase(end-1)]*fs/(2*pi);
frequency = movingAverage(frequency,round(smoothSec*fs));
end

function value = dualResidual(component)
[~,~,nHarmonic] = size(component);
residual = zeros(nHarmonic,1);
for ih = 1:nHarmonic
    difference = angle(component(:,2,ih).*conj(component(:,1,ih)));
    center = angle(mean(exp(1i*difference)));
    residual(ih) = sqrt(mean(angle(exp(1i*(difference-center))).^2));
end
value = rad2deg(mean(residual));
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
