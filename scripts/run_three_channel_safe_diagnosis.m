function diagnosis = run_three_channel_safe_diagnosis(healthyResult,faultResult)
%RUN_THREE_CHANNEL_SAFE_DIAGNOSIS BL600/PF50 feature comparison.
% Uses only the encoder-hidden common phase and the cross-validated safe
% path output produced by run_real_three_channel_phase.  The theoretical
% planet-tooth double-contact frequency is 2*f_mesh/Zp.

setup_paths;
if nargin < 1 || strlength(string(healthyResult))==0
    healthyResult = fullfile('results', ...
        'three_channel_BL_600RPM_Test01_26_05_09_14_27_04_MAT.mat');
end
if nargin < 2 || strlength(string(faultResult))==0
    faultResult = fullfile('results', ...
        'three_channel_PF50_600RPM_Test01_26_05_09_17_01_28_MAT.mat');
end

healthy = loadReport(healthyResult);
fault = loadReport(faultResult);
record = {healthy;fault};
conditionName = ["BL600_healthy";"PF50_planet_fault"];
methodName = ["single_selected_sensor";"naive_selected_pair"; ...
    "safe_aligned_common";"safe_aligned_differential"; ...
    "energy_preserving_common_diff"];
methodDisplay = ["Single sensor";"Naive pair";"Aligned common"; ...
    "Aligned differential";"Common+diff energy"];
nCondition = 2;
nMethod = numel(methodName);
planetToothCount = 31;

frequency = cell(nCondition,1);
spectrum = cell(nCondition,nMethod);
envelope = cell(nCondition,nMethod);
faultFrequencyHz = zeros(nCondition,1);
singleContactFrequencyHz = zeros(nCondition,1);
lineAmplitude = zeros(nCondition,nMethod,4);
lineToBackgroundDb = zeros(nCondition,nMethod);
rank1EnergyFraction = zeros(nCondition,nMethod);
rank3EnergyFraction = zeros(nCondition,nMethod);
faultCycleRepeatability = zeros(nCondition,nMethod);
faultCycleMatrix = cell(nCondition,nMethod);
singularValue = cell(nCondition,nMethod);
selectedSensorSubset = strings(nCondition,1);

for ic = 1:nCondition
    report = record{ic};
    [envelope(ic,:),phaseUse,selectedSensorSubset(ic)] = ...
        buildMethodEnvelopes(report);
    faultFrequencyHz(ic) = 2*report.blindMeshFrequencyHz/ ...
        planetToothCount;
    singleContactFrequencyHz(ic) = report.blindMeshFrequencyHz/ ...
        planetToothCount;
    for im = 1:nMethod
        [frequency{ic},spectrum{ic,im}] = oneSidedSpectrum( ...
            envelope{ic,im},report.fs);
        for ih = 1:4
            target = ih*faultFrequencyHz(ic);
            useLine = abs(frequency{ic}-target)<=0.35;
            lineAmplitude(ic,im,ih) = max(spectrum{ic,im}(useLine));
        end
        useBackground = frequency{ic}>=2 & frequency{ic}<=60;
        for ih = 1:4
            useBackground = useBackground & abs(frequency{ic}- ...
                ih*faultFrequencyHz(ic))>0.9;
        end
        lineRms = sqrt(mean(lineAmplitude(ic,im,:).^2,3));
        lineToBackgroundDb(ic,im) = 20*log10(lineRms/ ...
            median(spectrum{ic,im}(useBackground)));

        faultPhase = 2*phaseUse/planetToothCount;
        phaseMatrix = pg_phase_synchronous_matrix( ...
            envelope{ic,im},faultPhase,512);
        matrix = phaseMatrix.matrix(:,:,1);
        matrix = matrix-mean(matrix,2);
        faultCycleMatrix{ic,im} = matrix;
        [~,s,~] = svd(matrix,'econ');
        sv = diag(s);
        singularValue{ic,im} = sv;
        rank1EnergyFraction(ic,im) = sv(1)^2/sum(sv.^2);
        rank3EnergyFraction(ic,im) = ...
            sum(sv(1:min(3,numel(sv))).^2)/sum(sv.^2);
        faultCycleRepeatability(ic,im) = ...
            meanOffDiagonalCorrelation(matrix);
    end
end

lineSnrContrastDb = lineToBackgroundDb(2,:)-lineToBackgroundDb(1,:);
repeatabilityContrast = faultCycleRepeatability(2,:)- ...
    faultCycleRepeatability(1,:);
conditionColumn = strings(nCondition*nMethod,1);
methodColumn = strings(nCondition*nMethod,1);
sensorSubsetColumn = strings(nCondition*nMethod,1);
faultFrequencyColumn = zeros(nCondition*nMethod,1);
line1 = zeros(nCondition*nMethod,1);
line2 = zeros(nCondition*nMethod,1);
line3 = zeros(nCondition*nMethod,1);
line4 = zeros(nCondition*nMethod,1);
lineSnr = zeros(nCondition*nMethod,1);
rank1 = zeros(nCondition*nMethod,1);
rank3 = zeros(nCondition*nMethod,1);
repeatability = zeros(nCondition*nMethod,1);
row = 0;
for ic = 1:nCondition
    for im = 1:nMethod
        row = row+1;
        conditionColumn(row) = conditionName(ic);
        methodColumn(row) = methodName(im);
        sensorSubsetColumn(row) = selectedSensorSubset(ic);
        faultFrequencyColumn(row) = faultFrequencyHz(ic);
        line1(row) = lineAmplitude(ic,im,1);
        line2(row) = lineAmplitude(ic,im,2);
        line3(row) = lineAmplitude(ic,im,3);
        line4(row) = lineAmplitude(ic,im,4);
        lineSnr(row) = lineToBackgroundDb(ic,im);
        rank1(row) = rank1EnergyFraction(ic,im);
        rank3(row) = rank3EnergyFraction(ic,im);
        repeatability(row) = faultCycleRepeatability(ic,im);
    end
end
metrics = table(conditionColumn,methodColumn,sensorSubsetColumn, ...
    faultFrequencyColumn,line1,line2,line3,line4,lineSnr,rank1,rank3, ...
    repeatability,'VariableNames',{'condition','method','sensor_subset', ...
    'planet_fault_frequency_hz','fault_1x_amplitude', ...
    'fault_2x_amplitude','fault_3x_amplitude','fault_4x_amplitude', ...
    'fault_family_to_background_db','tpsvd_rank1_energy_fraction', ...
    'tpsvd_rank3_energy_fraction','fault_cycle_repeatability'});
contrast = table(methodName,lineSnrContrastDb(:), ...
    repeatabilityContrast(:),'VariableNames',{'method', ...
    'PF50_minus_BL_line_snr_db', ...
    'PF50_minus_BL_cycle_repeatability'});
writetable(metrics,fullfile('results', ...
    'three_channel_safe_diagnosis_metrics.csv'));
writetable(contrast,fullfile('results', ...
    'three_channel_safe_diagnosis_contrast.csv'));

color = [0.15 0.15 0.15;0.20 0.50 0.82;0.86 0.22 0.12; ...
    0.52 0.30 0.72;0.15 0.62 0.40];
fig = figure('Visible','off','Color','w','Position',[40 35 1240 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for ic = 1:nCondition
    nexttile;
    for im = 1:nMethod
        plot(frequency{ic},spectrum{ic,im},'Color',color(im,:), ...
            'LineWidth',0.9); hold on;
    end
    for ih = 1:4
        xline(ih*faultFrequencyHz(ic),'k:',sprintf('%d f_{pf}',ih), ...
            'HandleVisibility','off');
    end
    xline(singleContactFrequencyHz(ic),'Color',[0.55 0.55 0.55], ...
        'LineStyle','--','Label','f_m/Z_p','HandleVisibility','off');
    xlim([0 60]); grid on; box off;
    xlabel('Envelope frequency (Hz)'); ylabel('Amplitude');
    title(sprintf('%s, selected sensors %s',conditionName(ic), ...
        selectedSensorSubset(ic)),'Interpreter','none');
    legend(methodDisplay,'Location','best');
end
nexttile;
methodCategory = categorical(methodDisplay,methodDisplay,'Ordinal',true);
bar(methodCategory, ...
    lineToBackgroundDb.'); grid on; box off;
ylabel('Fault-family/background (dB)');
title('Theoretical planet-fault family');
legend(conditionName,'Interpreter','none','Location','best');
nexttile;
yyaxis left;
bar(methodCategory, ...
    lineSnrContrastDb,'FaceColor',[0.30 0.60 0.78]);
ylabel('PF50 - BL line SNR (dB)');
yyaxis right;
plot(methodCategory, ...
    repeatabilityContrast,'ro-','LineWidth',1.2,'MarkerFaceColor','r');
ylabel('PF50 - BL repeatability');
grid on; box off; title('Fault/healthy feature contrast');
exportgraphics(fig,fullfile('results', ...
    'three_channel_safe_diagnosis_spectrum.png'),'Resolution',190);
close(fig);

fig = figure('Visible','off','Color','w','Position',[55 45 1240 820]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
imagesc(linspace(0,360,size(faultCycleMatrix{2,2},2)), ...
    1:size(faultCycleMatrix{2,2},1),faultCycleMatrix{2,2}); axis xy;
xlabel('Planet-fault phase (deg)'); ylabel('Cycle index');
title('PF50 naive selected-pair cycles'); colorbar;
nexttile;
imagesc(linspace(0,360,size(faultCycleMatrix{2,5},2)), ...
    1:size(faultCycleMatrix{2,5},1),faultCycleMatrix{2,5}); axis xy;
xlabel('Planet-fault phase (deg)'); ylabel('Cycle index');
title('PF50 energy-preserving common+diff cycles'); colorbar;
nexttile;
for im = 1:nMethod
    cumulative = cumsum(singularValue{1,im}.^2)/ ...
        sum(singularValue{1,im}.^2);
    plot(1:numel(cumulative),cumulative,'Color',color(im,:), ...
        'LineWidth',1.0); hold on;
end
xlim([1 15]); ylim([0 1.02]); grid on; box off;
xlabel('Singular index'); ylabel('Cumulative energy');
title('BL600 fault-phase TPSVD');
legend(methodDisplay,'Location','best');
nexttile;
for im = 1:nMethod
    cumulative = cumsum(singularValue{2,im}.^2)/ ...
        sum(singularValue{2,im}.^2);
    plot(1:numel(cumulative),cumulative,'Color',color(im,:), ...
        'LineWidth',1.0); hold on;
end
xlim([1 15]); ylim([0 1.02]); grid on; box off;
xlabel('Singular index'); ylabel('Cumulative energy');
title('PF50 fault-phase TPSVD');
legend(methodDisplay,'Location','best');
exportgraphics(fig,fullfile('results', ...
    'three_channel_safe_diagnosis_tpsvd.png'),'Resolution',190);
close(fig);

diagnosis.healthyResult = healthyResult;
diagnosis.faultResult = faultResult;
diagnosis.condition = conditionName;
diagnosis.method = methodName;
diagnosis.selectedSensorSubset = selectedSensorSubset;
diagnosis.faultFrequencyHz = faultFrequencyHz;
diagnosis.singleContactFrequencyHz = singleContactFrequencyHz;
diagnosis.metrics = metrics;
diagnosis.contrast = contrast;
diagnosis.envelope = envelope;
diagnosis.frequency = frequency;
diagnosis.spectrum = spectrum;
diagnosis.faultCycleMatrix = faultCycleMatrix;
diagnosis.singularValue = singularValue;
diagnosis.interpretation = ['Known-condition feature comparison only. ', ...
    'The path topology is selected without labels or encoder samples; ', ...
    'classification accuracy requires additional repeated records.'];
save(fullfile('results','three_channel_safe_diagnosis.mat'), ...
    'diagnosis','-v7.3');

fprintf('Three-channel safe diagnosis complete.\n');
fprintf('  BL/PF50 planet-fault frequencies: %.4f / %.4f Hz\n', ...
    faultFrequencyHz(1),faultFrequencyHz(2));
disp(metrics);
disp(contrast);
end

function report = loadReport(file)
loaded = load(file,'report');
report = loaded.report;
required = {'safePath','adaptiveThree','rawSignal','harmonicOrder'};
for i = 1:numel(required)
    if ~isfield(report,required{i})
        error('run_three_channel_safe_diagnosis:ReportField', ...
            'Missing report.%s in %s.',required{i},file);
    end
end
end

function [envelopeSet,meshPhase,sensorSubset] = ...
    buildMethodEnvelopes(report)
sensorIndex = report.safePath.selectedSensorIndex;
sensorSubset = strjoin(string(sensorIndex),'+');
centerFrequencyHz = report.harmonicOrder*report.blindMeshFrequencyHz;
halfBandwidthHz = 30+4*(report.harmonicOrder-1);
raw = pg_extract_analytic_harmonics( ...
    report.rawSignal(:,sensorIndex),report.fs, ...
    centerFrequencyHz,halfBandwidthHz);
aligned = pg_extract_analytic_harmonics( ...
    report.safePath.correctedSignal,report.fs, ...
    centerFrequencyHz,halfBandwidthHz);
edge = round(0.50*report.fs);
use = (1+edge):(size(report.rawSignal,1)-edge);
rawComponent = raw.component(use,:,:);
alignedComponent = aligned.component(use,:,:);
meshPhase = report.adaptiveThree.meshPhase(use);
nSample = numel(use);
nHarmonic = numel(report.harmonicOrder);
singleEnvelope = zeros(nSample,1);
naiveEnvelope = zeros(nSample,1);
alignedEnvelope = zeros(nSample,1);
differentialEnvelope = zeros(nSample,1);
energyPreservingEnvelope = zeros(nSample,1);
for ih = 1:nHarmonic
    scale = sqrt(mean(abs(rawComponent(:,:,ih)).^2,'all'));
    singleEnvelope = singleEnvelope+ ...
        abs(rawComponent(:,1,ih))/max(scale,eps);
    naiveEnvelope = naiveEnvelope+abs(mean( ...
        rawComponent(:,:,ih),2))/max(scale,eps);
    alignedCommon = mean(alignedComponent(:,:,ih),2);
    alignedResidual = alignedComponent(:,:,ih)-alignedCommon;
    residualRms = sqrt(mean(abs(alignedResidual).^2,2));
    alignedEnvelope = alignedEnvelope+abs(alignedCommon)/max(scale,eps);
    differentialEnvelope = differentialEnvelope+ ...
        residualRms/max(scale,eps);
    % The common/differential transform is energy preserving.  This output
    % cannot erase a localized channel feature when the coherent common
    % mode alone is weak.
    energyPreservingEnvelope = energyPreservingEnvelope+sqrt( ...
        abs(alignedCommon).^2+residualRms.^2)/max(scale,eps);
end
singleEnvelope = detrend(singleEnvelope/nHarmonic);
naiveEnvelope = detrend(naiveEnvelope/nHarmonic);
alignedEnvelope = detrend(alignedEnvelope/nHarmonic);
differentialEnvelope = detrend(differentialEnvelope/nHarmonic);
energyPreservingEnvelope = detrend( ...
    energyPreservingEnvelope/nHarmonic);
envelopeSet = {singleEnvelope,naiveEnvelope,alignedEnvelope, ...
    differentialEnvelope,energyPreservingEnvelope};
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
correlation = corrcoef(matrix.');
use = triu(true(size(correlation)),1);
value = mean(correlation(use),'omitnan');
end
