function diagnosis = run_real_phase_aligned_diagnosis
%RUN_REAL_PHASE_ALIGNED_DIAGNOSIS Fault-feature comparison after alignment.
% This script uses the frozen encoder-hidden phase result. It diagnoses the
% known 50% sun-tooth-break record; it does not estimate classification
% accuracy because no healthy/25% records are currently available.

setup_paths;
loaded = load(fullfile('results','real_tacholess_joint_phase.mat'),'report');
phaseReport = loaded.report;
fs = phaseReport.fs;
t = phaseReport.t;
meshFrequencyBlindHz = phaseReport.blindMeshFrequencyHz;
sunToothCount = 21;
faultFrequencyBlindHz = 3*meshFrequencyBlindHz/sunToothCount;
faultFrequencyEncoderHz = 3*mean(phaseReport.truth.meshFrequencyHz)/ ...
    sunToothCount;

alignedComponent = phaseReport.hybrid.correctedComponent;
% corrected = measured*conj(path), so this exactly recovers the measured
% analytic components without rereading or filtering the TDMS file.
measuredComponent = alignedComponent.*phaseReport.hybrid.pathPhasor;
[nSample,~,nHarmonic] = size(measuredComponent);

singleEnvelope = zeros(nSample,1);
naiveEnvelope = zeros(nSample,1);
alignedEnvelope = zeros(nSample,1);
for ih = 1:nHarmonic
    z0 = measuredComponent(:,1,ih);
    z90 = measuredComponent(:,2,ih);
    zSingle = z0;
    zNaive = (z0+z90)/2;
    zAligned = (alignedComponent(:,1,ih)+ ...
        alignedComponent(:,2,ih))/2;
    commonScale = sqrt(mean((abs(z0).^2+abs(z90).^2)/2));
    singleEnvelope = singleEnvelope+abs(zSingle)/commonScale;
    naiveEnvelope = naiveEnvelope+abs(zNaive)/commonScale;
    alignedEnvelope = alignedEnvelope+abs(zAligned)/commonScale;
end
singleEnvelope = detrend(singleEnvelope/nHarmonic);
naiveEnvelope = detrend(naiveEnvelope/nHarmonic);
alignedEnvelope = detrend(alignedEnvelope/nHarmonic);

method = {'single_0deg';'naive_dual_fusion';'phase_aligned_fusion'};
envelopeSet = {singleEnvelope,naiveEnvelope,alignedEnvelope};
nMethod = numel(method);
frequency = [];
spectrumSet = cell(nMethod,1);
faultHarmonicAmplitude = zeros(nMethod,4);
faultToBackgroundDb = zeros(nMethod,1);
rank1EnergyFraction = zeros(nMethod,1);
rank3EnergyFraction = zeros(nMethod,1);
carrierOrder12Amplitude = zeros(nMethod,1);
cycleRepeatability = zeros(nMethod,1);
cycleMatrix = cell(nMethod,1);
tpsvdPattern = cell(nMethod,1);
singularValue = cell(nMethod,1);
orderSpectrum = cell(nMethod,1);

for im = 1:nMethod
    [frequency,spectrumSet{im}] = oneSidedSpectrum(envelopeSet{im},fs);
    for ihf = 1:4
        target = ihf*faultFrequencyBlindHz;
        useLine = abs(frequency-target)<=0.55;
        faultHarmonicAmplitude(im,ihf) = max(spectrumSet{im}(useLine));
    end
    useBackground = frequency>=5 & frequency<=110;
    for ihf = 1:4
        useBackground = useBackground & ...
            abs(frequency-ihf*faultFrequencyBlindHz)>1.2;
    end
    lineRms = sqrt(mean(faultHarmonicAmplitude(im,:).^2));
    faultToBackgroundDb(im) = 20*log10(lineRms/ ...
        median(spectrumSet{im}(useBackground)));

    phaseMatrix = pg_phase_synchronous_matrix(envelopeSet{im}, ...
        phaseReport.hybrid.carrierPhase,2048);
    matrix = phaseMatrix.matrix(:,:,1);
    matrix = matrix-mean(matrix,2);
    cycleMatrix{im} = matrix;
    [u,s,v] = svd(matrix,'econ');
    sv = diag(s);
    singularValue{im} = sv;
    rank1EnergyFraction(im) = sv(1)^2/sum(sv.^2);
    rank3EnergyFraction(im) = sum(sv(1:min(3,numel(sv))).^2)/sum(sv.^2);
    rankUse = 1:min(3,numel(sv));
    reconstruction = u(:,rankUse)*s(rankUse,rankUse)*v(:,rankUse)';
    pattern = mean(reconstruction,1).';
    tpsvdPattern{im} = pattern;
    patternSpectrum = 2*abs(fft(pattern.*hann(numel(pattern))))/ ...
        sum(hann(numel(pattern)));
    orderSpectrum{im} = patternSpectrum(1:floor(numel(pattern)/2));
    [~,order12Index] = min(abs((0:numel(orderSpectrum{im})-1)'-12));
    carrierOrder12Amplitude(im) = orderSpectrum{im}(order12Index);
    cycleRepeatability(im) = meanOffDiagonalCorrelation(matrix);
end

faultGainVsNaiveDb = 20*log10( ...
    sqrt(mean(faultHarmonicAmplitude(:,1:4).^2,2))/ ...
    sqrt(mean(faultHarmonicAmplitude(2,1:4).^2)));
metrics = table(method,faultHarmonicAmplitude(:,1), ...
    faultHarmonicAmplitude(:,2),faultHarmonicAmplitude(:,3), ...
    faultHarmonicAmplitude(:,4),faultToBackgroundDb,faultGainVsNaiveDb, ...
    rank1EnergyFraction,rank3EnergyFraction,carrierOrder12Amplitude, ...
    cycleRepeatability, ...
    'VariableNames',{'method','fault_1x_amplitude','fault_2x_amplitude', ...
    'fault_3x_amplitude','fault_4x_amplitude','fault_to_background_db', ...
    'fault_gain_vs_naive_db','tpsvd_rank1_energy_fraction', ...
    'tpsvd_rank3_energy_fraction','carrier_order12_amplitude', ...
    'cycle_repeatability'});
writetable(metrics,fullfile('results', ...
    'real_phase_aligned_diagnosis_metrics.csv'));

fig = figure('Visible','off','Color','w','Position',[50 40 1220 830]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
color = [0.25 0.25 0.25;0.20 0.45 0.80;0.85 0.20 0.10];
nexttile;
for im = 1:nMethod
    plot(t,envelopeSet{im},'Color',color(im,:),'LineWidth',0.8); hold on;
end
xlim([0 min(1.5,t(end))]); grid on; box off;
xlabel('Time (s)'); ylabel('Normalized multiharmonic envelope');
title('Measured 50% sun-tooth break');
legend(method,'Interpreter','none','Location','best');

nexttile;
for im = 1:nMethod
    plot(frequency,spectrumSet{im},'Color',color(im,:),'LineWidth',1.0); hold on;
end
for ihf = 1:4
    xline(ihf*faultFrequencyEncoderHz,'k:',sprintf('%d f_{sf}',ihf), ...
        'HandleVisibility','off');
end
xlim([0 105]); grid on; box off;
xlabel('Envelope frequency (Hz)'); ylabel('Amplitude');
title(sprintf('Sun-fault family, f_{sf}=%.2f Hz',faultFrequencyEncoderHz));
legend(method,'Interpreter','none','Location','best');

nexttile;
bar(1:4,faultHarmonicAmplitude.'); grid on; box off;
xticks(1:4); xticklabels({'1f_{sf}','2f_{sf}','3f_{sf}','4f_{sf}'});
xlabel('Fault harmonic'); ylabel('Envelope amplitude');
title('Fault-family amplitudes');
legend(method,'Interpreter','none','Location','best');

nexttile;
yyaxis left;
bar(categorical({'Single','Naive dual','Aligned dual'}, ...
    {'Single','Naive dual','Aligned dual'},'Ordinal',true), ...
    faultToBackgroundDb,'FaceColor',[0.25 0.55 0.80]);
ylabel('Fault/background ratio (dB)');
yyaxis right;
plot(categorical({'Single','Naive dual','Aligned dual'}, ...
    {'Single','Naive dual','Aligned dual'},'Ordinal',true), ...
    cycleRepeatability,'ro-','LineWidth',1.2,'MarkerFaceColor','r');
ylabel('Carrier-cycle repeatability'); ylim([0 1]);
grid on; box off; title('Diagnostic visibility and repeatability');
exportgraphics(fig,fullfile('results', ...
    'real_phase_aligned_diagnosis_spectrum.png'),'Resolution',190);
close(fig);

fig = figure('Visible','off','Color','w','Position',[60 50 1220 830]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
imagesc(linspace(0,360,size(cycleMatrix{2},2)), ...
    1:size(cycleMatrix{2},1),cycleMatrix{2}); axis xy;
xlabel('Carrier angle (deg)'); ylabel('Cycle index');
title('Naive dual-fusion carrier cycles'); colorbar;
nexttile;
imagesc(linspace(0,360,size(cycleMatrix{3},2)), ...
    1:size(cycleMatrix{3},1),cycleMatrix{3}); axis xy;
xlabel('Carrier angle (deg)'); ylabel('Cycle index');
title('Phase-aligned carrier cycles'); colorbar;

nexttile;
for im = 1:nMethod
    cumulative = cumsum(singularValue{im}.^2)/sum(singularValue{im}.^2);
    plot(1:numel(cumulative),cumulative,'o-','Color',color(im,:), ...
        'LineWidth',1.0,'MarkerSize',3); hold on;
end
xlim([1 min(12,numel(singularValue{1}))]); ylim([0 1.02]);
grid on; box off; xlabel('Singular index');
ylabel('Cumulative energy'); title('Carrier-cycle TPSVD');
legend(method,'Interpreter','none','Location','best');

nexttile;
carrierOrder = (0:numel(orderSpectrum{1})-1).';
for im = 1:nMethod
    plot(carrierOrder,orderSpectrum{im},'Color',color(im,:), ...
        'LineWidth',1.0); hold on;
end
xline(12,'k:','12 = f_{sf}/f_c','HandleVisibility','off');
xline(24,'k:','24','HandleVisibility','off');
xlim([0 40]); grid on; box off;
xlabel('Carrier order'); ylabel('Rank-3 pattern amplitude');
title('TPSVD fault order spectrum');
legend(method,'Interpreter','none','Location','best');
exportgraphics(fig,fullfile('results', ...
    'real_phase_aligned_tpsvd.png'),'Resolution',190);
close(fig);

diagnosis.sourcePhaseResult = fullfile('results', ...
    'real_tacholess_joint_phase.mat');
diagnosis.t = t;
diagnosis.faultFrequencyBlindHz = faultFrequencyBlindHz;
diagnosis.faultFrequencyEncoderHz = faultFrequencyEncoderHz;
diagnosis.metrics = metrics;
diagnosis.envelope = [singleEnvelope,naiveEnvelope,alignedEnvelope];
diagnosis.frequency = frequency;
diagnosis.spectrum = [spectrumSet{1},spectrumSet{2},spectrumSet{3}];
diagnosis.cycleMatrix = cycleMatrix;
diagnosis.tpsvdPattern = tpsvdPattern;
diagnosis.singularValue = singularValue;
diagnosis.interpretation = ['Feature visibility for one known 50% sun ', ...
    'tooth-break record; not a classification-accuracy experiment.'];
save(fullfile('results','real_phase_aligned_diagnosis.mat'), ...
    'diagnosis','-v7.3');

fprintf('Real phase-aligned diagnosis complete.\n');
fprintf('  sun fault frequency blind/encoder %.4f / %.4f Hz\n', ...
    faultFrequencyBlindHz,faultFrequencyEncoderHz);
disp(metrics);
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
