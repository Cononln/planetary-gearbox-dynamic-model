function output = simulate_loaded_v2_responses
%SIMULATE_LOADED_V2_RESPONSES Generate provisional loaded-model figures.

setup_paths;
p = pg_parameters_loaded_v2;
model = pg_ring_modal_model(p);
fs = 51200;
discardTime = 0.20;
recordDuration = 3.00;
duration = discardTime+recordDuration;

cases = {pg_fault_case('none',0,p), ...
         pg_fault_case('sun',0.25,p), ...
         pg_fault_case('sun',0.50,p)};
caseNames = {'Healthy','Sun broken tooth 25%','Sun broken tooth 50%'};
colors = [0.05 0.05 0.05;0.05 0.35 0.85;0.85 0.15 0.10];
nCase = numel(cases);
responses = cell(1,nCase);

for ic = 1:nCase
    fprintf('Loaded-v2 simulation: %s\n',caseNames{ic});
    fullResponse = pg_simulate_loaded_v2_response(duration,p,cases{ic}, ...
        model,'fs',fs);
    keep = fullResponse.time>=discardTime;
    fullResponse.time = fullResponse.time(keep)-discardTime;
    fullResponse.acceleration = fullResponse.acceleration(keep,:);
    fullResponse.acceleration = fullResponse.acceleration- ...
        mean(fullResponse.acceleration,1);
    fullResponse.meshForceSunPlanet = ...
        fullResponse.meshForceSunPlanet(keep,:);
    fullResponse.meshForceRingPlanet = ...
        fullResponse.meshForceRingPlanet(keep,:);
    fullResponse.inputTorque = fullResponse.inputTorque(keep,:);
    assert(all(isfinite(fullResponse.acceleration),'all'));
    responses{ic} = fullResponse;
end

time = responses{1}.time;
n = numel(time);
rmsAcceleration = zeros(nCase,2);
crestFactor = zeros(nCase,2);
kurtosisValue = zeros(nCase,2);
meanCompressionSp = zeros(nCase,p.model.nPlanet);
for ic = 1:nCase
    x = responses{ic}.acceleration;
    rmsAcceleration(ic,:) = sqrt(mean(x.^2,1));
    crestFactor(ic,:) = max(abs(x),[],1)./rmsAcceleration(ic,:);
    variance = mean(x.^2,1);
    kurtosisValue(ic,:) = mean(x.^4,1)./(variance.^2);
    meanCompressionSp(ic,:) = -mean(responses{ic}.meshForceSunPlanet,1);
end

window = 0.5-0.5*cos(2*pi*(0:n-1)'/(n-1));
frequency = (0:floor(n/2))'*fs/n;
spectrum = zeros(numel(frequency),nCase,2);
for ic = 1:nCase
    for is = 1:2
        X = fft(responses{ic}.acceleration(:,is).*window);
        spectrum(:,ic,is) = 2*abs(X(1:numel(frequency)))/sum(window);
    end
end

output.parameters = p;
output.ringModel = model;
output.caseNames = caseNames;
output.responses = responses;
output.frequency = frequency;
output.spectrum = spectrum;
output.rmsAcceleration = rmsAcceleration;
output.crestFactor = crestFactor;
output.kurtosis = kurtosisValue;
output.meanCompressionSunPlanet = meanCompressionSp;
save(fullfile('results','loaded_v2_response_cases.mat'),'output','-v7.3');

signalTable = table(time, ...
    responses{1}.acceleration(:,1),responses{1}.acceleration(:,2), ...
    responses{2}.acceleration(:,1),responses{2}.acceleration(:,2), ...
    responses{3}.acceleration(:,1),responses{3}.acceleration(:,2), ...
    'VariableNames',{'time_s','healthy_0deg_ms2','healthy_90deg_ms2', ...
    'sun25_0deg_ms2','sun25_90deg_ms2','sun50_0deg_ms2','sun50_90deg_ms2'});
writetable(signalTable,fullfile('results','loaded_v2_signals_3s.csv'));

metricsTable = table(string(caseNames(:)), ...
    rmsAcceleration(:,1),rmsAcceleration(:,2), ...
    crestFactor(:,1),crestFactor(:,2), ...
    kurtosisValue(:,1),kurtosisValue(:,2), ...
    'VariableNames',{'case_name','rms_0deg_ms2','rms_90deg_ms2', ...
    'crest_0deg','crest_90deg','kurtosis_0deg','kurtosis_90deg'});
writetable(metricsTable,fullfile('results','loaded_v2_metrics.csv'));

% One-second raw waveforms.
showTime = time<=1;
fig = figure('Visible','off','Color','w','Position',[60 40 1180 850]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
for ic = 1:nCase
    for is = 1:2
        nexttile;
        plot(time(showTime),responses{ic}.acceleration(showTime,is), ...
            'Color',colors(ic,:),'LineWidth',0.65);
        xlim([0 1]); grid on; box off;
        xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
        title(sprintf('%s, %d deg: RMS %.2f, crest %.2f', ...
            caseNames{ic},responses{ic}.sensorAnglesDeg(is), ...
            rmsAcceleration(ic,is),crestFactor(ic,is)));
    end
end
sgtitle(sprintf(['Loaded 18+6 DOF model, carrier load %.0f N m ' ...
    '(provisional)'],p.operating.outputTorque));
exportgraphics(fig,fullfile('results','loaded_v2_time_1s.png'), ...
    'Resolution',190);
close(fig);

% Mean sun-planet contact force: a loaded signal instead of zero-mean TE force.
forceTime = time<=0.50;
fig = figure('Visible','off','Color','w','Position',[80 80 1120 660]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
for panel = 1:2
    nexttile;
    for ic = 1:nCase
        compression = -responses{ic}.meshForceSunPlanet(:,1);
        if panel == 1
            plot(time(forceTime),compression(forceTime), ...
                'Color',colors(ic,:),'LineWidth',0.8); hold on;
        else
            dynamicCompression = compression-mean(compression);
            plot(time(forceTime),dynamicCompression(forceTime), ...
                'Color',colors(ic,:),'LineWidth',0.8); hold on;
        end
    end
    xlim([0 0.5]); grid on; box off;
    xlabel('Time (s)');
    if panel == 1
        ylabel('Compression force (N)');
        title('Planet 1 sun-mesh loaded contact force');
    else
        ylabel('Dynamic force (N)');
        title('Mean-removed contact force');
    end
    legend(caseNames,'Location','best');
end
exportgraphics(fig,fullfile('results','loaded_v2_mesh_force_0p5s.png'), ...
    'Resolution',190);
close(fig);

% Main low-frequency spectrum: linear physical amplitude starting at zero.
inLow = frequency<=500;
fig = figure('Visible','off','Color','w','Position',[80 80 1120 680]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
for is = 1:2
    nexttile;
    for ic = 1:nCase
        plot(frequency(inLow),spectrum(inLow,ic,is), ...
            'Color',colors(ic,:),'LineWidth',1.0); hold on;
    end
    xline(p.kin.f_mesh,'--','f_m=168 Hz','Color',[0.3 0.3 0.3]);
    lowMaximum = max(spectrum(inLow,:,is),[],'all');
    xlim([0 500]); ylim([0 1.05*lowMaximum]); grid on; box off;
    ax = gca; ax.YAxis.Exponent = 0; ytickformat('%.4g');
    xlabel('Frequency (Hz)'); ylabel('Amplitude (m/s^2)');
    title(sprintf('Sensor %d deg, Delta f = %.3f Hz', ...
        responses{1}.sensorAnglesDeg(is),frequency(2)-frequency(1)));
    legend(caseNames,'Location','northeast');
end
exportgraphics(fig,fullfile('results','loaded_v2_spectrum_0_500Hz.png'), ...
    'Resolution',190);
close(fig);

% Mesh-band detail: linear physical amplitude starting at zero.
inMesh = frequency>=140 & frequency<=200;
fig = figure('Visible','off','Color','w','Position',[80 80 1120 680]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
for is = 1:2
    nexttile;
    for ic = 1:nCase
        plot(frequency(inMesh),spectrum(inMesh,ic,is), ...
            'Color',colors(ic,:), ...
            'LineWidth',1.0); hold on;
    end
    for order = -3:3
        xline(p.kin.f_mesh+order*p.kin.f_c,':', ...
            'Color',[0.65 0.65 0.65],'HandleVisibility','off');
    end
    meshMaximum = max(spectrum(inMesh,:,is),[],'all');
    xlim([140 200]); ylim([0 1.05*meshMaximum]); grid on; box off;
    ax = gca; ax.YAxis.Exponent = 0; ytickformat('%.4g');
    xlabel('Frequency (Hz)'); ylabel('Amplitude (m/s^2)');
    title(sprintf('Sensor %d deg: mesh band and carrier orders', ...
        responses{1}.sensorAnglesDeg(is)));
    legend(caseNames,'Location','southwest');
end
exportgraphics(fig,fullfile('results','loaded_v2_mesh_band_140_200Hz.png'), ...
    'Resolution',190);
close(fig);

% Paper-style full spectrum: one case per panel, all amplitude axes start at 0.
inFull = frequency<=5500;
fullMaximum = zeros(1,2);
for is = 1:2
    fullMaximum(is) = max(spectrum(inFull,:,is),[],'all');
end
fig = figure('Visible','off','Color','w','Position',[60 40 1180 850]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
for ic = 1:nCase
    for is = 1:2
        nexttile;
        plot(frequency(inFull),spectrum(inFull,ic,is), ...
            'Color',colors(ic,:),'LineWidth',0.75);
        xlim([0 5500]); ylim([0 1.05*fullMaximum(is)]); grid on; box off;
        xlabel('Frequency (Hz)'); ylabel('Amplitude (m/s^2)');
        title(sprintf('%s, sensor %d deg',caseNames{ic}, ...
            responses{1}.sensorAnglesDeg(is)));
    end
end
exportgraphics(fig,fullfile('results','loaded_v2_full_spectrum_0_5500Hz.png'), ...
    'Resolution',190);
close(fig);

% Envelope spectrum around the first elastic-ring pair.
frequencyFull = (0:n-1)'*fs/n;
firstBand = [1700 2050];
bandMask = (frequencyFull>=firstBand(1) & frequencyFull<=firstBand(2)) | ...
    (frequencyFull>=fs-firstBand(2) & frequencyFull<=fs-firstBand(1));
envelopeSpectrum = zeros(numel(frequency),nCase,2);
for ic = 1:nCase
    for is = 1:2
        bandSignal = real(ifft(fft(responses{ic}.acceleration(:,is)).* ...
            bandMask));
        envelope = abs(hilbert(bandSignal));
        envelope = envelope-mean(envelope);
        E = fft(envelope.*window);
        envelopeSpectrum(:,ic,is) = ...
            2*abs(E(1:numel(frequency)))/sum(window);
    end
end
inEnvelope = frequency<=120;
fig = figure('Visible','off','Color','w','Position',[80 80 1120 680]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
for is = 1:2
    nexttile;
    for ic = 1:nCase
        plot(frequency(inEnvelope),envelopeSpectrum(inEnvelope,ic,is), ...
            'Color',colors(ic,:),'LineWidth',1.0); hold on;
    end
    for harmonic = 1:5
        xline(harmonic*p.kin.f_sun_fault,':', ...
            'Color',[0.65 0.65 0.65],'HandleVisibility','off');
    end
    envelopeMaximum = max(envelopeSpectrum(inEnvelope,:,is),[],'all');
    xlim([0 120]); ylim([0 1.05*envelopeMaximum]); grid on; box off;
    xlabel('Frequency (Hz)'); ylabel('Envelope amplitude (m/s^2)');
    title(sprintf('Sensor %d deg: %.2f-%.2f kHz envelope', ...
        responses{1}.sensorAnglesDeg(is),firstBand(1)/1000,firstBand(2)/1000));
    legend(caseNames,'Location','northeast');
end
exportgraphics(fig,fullfile('results','loaded_v2_envelope_0_120Hz.png'), ...
    'Resolution',190);
close(fig);

fprintf('Loaded-v2 figures complete.\n');
for ic = 1:nCase
    fprintf(['  %-24s RMS %.3f/%.3f m/s^2; crest %.2f/%.2f; ' ...
        'kurtosis %.2f/%.2f\n'],caseNames{ic}, ...
        rmsAcceleration(ic,1),rmsAcceleration(ic,2), ...
        crestFactor(ic,1),crestFactor(ic,2), ...
        kurtosisValue(ic,1),kurtosisValue(ic,2));
end
end
