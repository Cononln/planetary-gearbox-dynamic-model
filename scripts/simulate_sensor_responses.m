function output = simulate_sensor_responses
%SIMULATE_SENSOR_RESPONSES Generate default 18-DOF LPM dual-sensor signals.

setup_paths;
p = pg_parameters;
discardTime = 0.10;
recordDuration = 1.00;
duration = discardTime+recordDuration;
cases = {pg_fault_case('none',0,p), ...
         pg_fault_case('sun',0.25,p), ...
         pg_fault_case('sun',0.50,p)};
caseNames = {'Healthy','Sun fault 25%','Sun fault 50%'};

responses = cell(1,numel(cases));
for ic = 1:numel(cases)
    fprintf('Simulating %s ...\n',caseNames{ic});
    fullResponse = pg_simulate_lpm_response(duration,p,cases{ic});
    keep = fullResponse.time>=discardTime;
    fullResponse.time = fullResponse.time(keep)-discardTime;
    fullResponse.acceleration = fullResponse.acceleration(keep,:);
    fullResponse.acceleration = fullResponse.acceleration- ...
        mean(fullResponse.acceleration,1);
    assert(all(isfinite(fullResponse.acceleration),'all'));
    responses{ic} = fullResponse;
end

time = responses{1}.time;
nCase = numel(cases);
rmsAcceleration = zeros(nCase,2);
peakAcceleration = zeros(nCase,2);
crestFactor = zeros(nCase,2);
for ic = 1:nCase
    signal = responses{ic}.acceleration;
    rmsAcceleration(ic,:) = sqrt(mean(signal.^2,1));
    peakAcceleration(ic,:) = max(abs(signal),[],1);
    crestFactor(ic,:) = peakAcceleration(ic,:)./rmsAcceleration(ic,:);
end

output.caseNames = caseNames;
output.responses = responses;
output.rmsAcceleration = rmsAcceleration;
output.peakAcceleration = peakAcceleration;
output.crestFactor = crestFactor;
output.units = 'm/s^2';
output.modelType = responses{1}.modelType;
save(fullfile('results','sensor_response_sun_cases_1s.mat'),'output','-v7.3');

responseTable = table(time, ...
    responses{1}.acceleration(:,1),responses{1}.acceleration(:,2), ...
    responses{2}.acceleration(:,1),responses{2}.acceleration(:,2), ...
    responses{3}.acceleration(:,1),responses{3}.acceleration(:,2), ...
    'VariableNames',{'time_s','healthy_0deg_ms2','healthy_90deg_ms2', ...
    'sun25_0deg_ms2','sun25_90deg_ms2','sun50_0deg_ms2','sun50_90deg_ms2'});
writetable(responseTable,fullfile('results','sensor_response_sun_cases_1s.csv'));

% Time-domain figure in mm/s^2 avoids hidden scientific-notation exponents.
displayScale = 1e3;
fig = figure('Visible','off','Color','w','Position',[80 80 1120 800]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
for ic = 1:nCase
    for is = 1:2
        nexttile;
        plot(time,displayScale*responses{ic}.acceleration(:,is), ...
            'Color',[0.05 0.25 0.70], ...
            'LineWidth',0.7);
        xlim([time(1),time(end)]);
        xlabel('Time (s)'); ylabel('Acceleration (mm/s^2)');
        title(sprintf('%s, sensor %d deg; RMS %.2f mm/s^2, crest %.2f', ...
            caseNames{ic},responses{ic}.sensorAnglesDeg(is), ...
            displayScale*rmsAcceleration(ic,is),crestFactor(ic,is)));
        grid on;
    end
end
exportgraphics(fig,fullfile('results','sensor_response_time_1s.png'),'Resolution',180);
close(fig);

% Spectral comparison for both sensor positions.
n = numel(time);
window = 0.5-0.5*cos(2*pi*(0:n-1)'/(n-1));
frequency = (0:floor(n/2))'*p.model.fs/n;
spectrum = zeros(numel(frequency),nCase,2);
for ic = 1:nCase
    for is = 1:2
        x = responses{ic}.acceleration(:,is).*window;
        X = fft(x);
        spectrum(:,ic,is) = 2*abs(X(1:numel(frequency)))/sum(window);
    end
end
% Log-scale diagnostic spectrum.
fig = figure('Visible','off','Color','w','Position',[100 100 1040 620]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
colors = [0 0 0;0.05 0.35 0.85;0.85 0.15 0.10];
for is = 1:2
    nexttile;
    curveHandle = gobjects(1,nCase);
    for ic = 1:nCase
        curveHandle(ic) = semilogy(frequency,displayScale*spectrum(:,ic,is), ...
            'Color',colors(ic,:), ...
            'LineWidth',1.0); hold on;
    end
    xlim([0 5500]); xlabel('Frequency (Hz)');
    ylabel('Acceleration amplitude (mm/s^2)');
    title(sprintf('Sensor %d deg spectrum',responses{1}.sensorAnglesDeg(is)));
    legend(curveHandle,caseNames,'Location','best'); grid on;
    xline(p.sensor.pathNaturalFrequencyHz,':', ...
        'Color',[0.45 0.45 0.45],'HandleVisibility','off');
end
exportgraphics(fig,fullfile('results','sensor_response_spectrum_log_1s.png'), ...
    'Resolution',180);
close(fig);

% Paper-style linear spectrum: separate panels prevent the severe case from
% hiding the healthy response, and the horizontal range follows the reference.
maxPlotFrequency = 10000;
inPlotBand = frequency<=maxPlotFrequency;
commonMaximum = zeros(1,2);
for is = 1:2
    values = spectrum(inPlotBand,:,is);
    commonMaximum(is) = max(values,[],'all');
end
fig = figure('Visible','off','Color','w','Position',[80 50 1120 820]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
for ic = 1:nCase
    for is = 1:2
        nexttile;
        plot(frequency(inPlotBand),displayScale*spectrum(inPlotBand,ic,is), ...
            'Color',colors(ic,:),'LineWidth',0.8);
        xlim([0 maxPlotFrequency]);
        ylim([0 1.05*displayScale*commonMaximum(is)]);
        xlabel('Frequency (Hz)');
        ylabel('Single-sided amplitude (mm/s^2)');
        title(sprintf('%s, sensor %d deg',caseNames{ic}, ...
            responses{ic}.sensorAnglesDeg(is)));
        box off;
    end
end
exportgraphics(fig,fullfile('results','sensor_response_spectrum_1s.png'), ...
    'Resolution',180);
close(fig);

% Resonance-band envelope spectrum exposes the 24 Hz sun-fault modulation.
envelopeMaxFrequency = 120;
firstBand = [1700 2100];
envelopeSpectrum = zeros(numel(frequency),nCase,2);
frequencyFull = (0:n-1)'*p.model.fs/n;
bandMask = (frequencyFull>=firstBand(1) & frequencyFull<=firstBand(2)) | ...
    (frequencyFull>=p.model.fs-firstBand(2) & ...
     frequencyFull<=p.model.fs-firstBand(1));
for ic = 1:nCase
    for is = 1:2
        x = responses{ic}.acceleration(:,is);
        bandSignal = real(ifft(fft(x).*bandMask));
        envelope = abs(hilbert(bandSignal));
        envelope = envelope-mean(envelope);
        E = fft(envelope.*window);
        envelopeSpectrum(:,ic,is) = ...
            2*abs(E(1:numel(frequency)))/sum(window);
    end
end
inEnvelopeBand = frequency<=envelopeMaxFrequency;
fig = figure('Visible','off','Color','w','Position',[100 100 1040 620]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
for is = 1:2
    nexttile;
    curveHandle = gobjects(1,nCase);
    for ic = 1:nCase
        curveHandle(ic) = plot(frequency(inEnvelopeBand), ...
            displayScale*envelopeSpectrum(inEnvelopeBand,ic,is), ...
            'Color',colors(ic,:),'LineWidth',1.0); hold on;
    end
    for harmonic = 1:5
        xline(harmonic*p.kin.f_sun_fault,':','Color',[0.45 0.45 0.45], ...
            'HandleVisibility','off');
    end
    xlim([0 envelopeMaxFrequency]); xlabel('Frequency (Hz)');
    ylabel('Envelope amplitude (mm/s^2)');
    title(sprintf('Sensor %d deg, %.1f-%.1f kHz envelope spectrum', ...
        responses{1}.sensorAnglesDeg(is),firstBand(1)/1000,firstBand(2)/1000));
    legend(curveHandle,caseNames,'Location','best'); box off;
end
exportgraphics(fig,fullfile('results','sensor_response_envelope_spectrum_1s.png'), ...
    'Resolution',180);
close(fig);

fprintf('Sensor-response simulation PASSED\n');
for ic = 1:nCase
    fprintf('  %-13s RMS 0/90 deg: %.3f / %.3f m/s^2; crest: %.2f / %.2f\n', ...
        caseNames{ic},rmsAcceleration(ic,1),rmsAcceleration(ic,2), ...
        crestFactor(ic,1),crestFactor(ic,2));
end
end
