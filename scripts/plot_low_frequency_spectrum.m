function report = plot_low_frequency_spectrum
%PLOT_LOW_FREQUENCY_SPECTRUM Plot 0-200 Hz and distinguish two sideband sets.

setup_paths;
p = pg_parameters;
loaded = load(fullfile('results','sensor_response_sun_cases_1s.mat'),'output');
output = loaded.output;
responses = output.responses;
caseNames = output.caseNames;
nCase = numel(responses);
fs = responses{1}.fs;
n = numel(responses{1}.time);
window = 0.5-0.5*cos(2*pi*(0:n-1)'/(n-1));
frequency = (0:floor(n/2))'*fs/n;
spectrum = zeros(numel(frequency),nCase,2);
complexSpectrum = complex(zeros(numel(frequency),nCase,2));

for ic = 1:nCase
    for is = 1:2
        x = responses{ic}.acceleration(:,is);
        x = x-mean(x);
        X = fft(x.*window);
        complexSpectrum(:,ic,is) = 2*X(1:numel(frequency))/sum(window);
        spectrum(:,ic,is) = abs(complexSpectrum(:,ic,is));
    end
end

targetFrequency = [p.kin.f_sun_fault,2*p.kin.f_sun_fault, ...
    3*p.kin.f_sun_fault,4*p.kin.f_sun_fault,5*p.kin.f_sun_fault, ...
    p.kin.f_mesh-p.kin.f_sun_fault, ...
    p.kin.f_mesh-p.model.nPlanet*p.kin.f_c, ...
    p.kin.f_mesh-2*p.kin.f_c,p.kin.f_mesh-p.kin.f_c, ...
    p.kin.f_mesh,p.kin.f_mesh+p.kin.f_c, ...
    p.kin.f_mesh+2*p.kin.f_c, ...
    p.kin.f_mesh+p.model.nPlanet*p.kin.f_c, ...
    p.kin.f_mesh+p.kin.f_sun_fault];
targetNames = {'24','48','72','96','120','144','162','164','166', ...
    '168','170','172','174','192'};
targetAmplitude = zeros(numel(targetFrequency),nCase,2);
targetPhaseDifferenceDeg = zeros(numel(targetFrequency),nCase);
for it = 1:numel(targetFrequency)
    [~,idx] = min(abs(frequency-targetFrequency(it)));
    targetAmplitude(it,:,:) = spectrum(idx,:,:);
    targetPhaseDifferenceDeg(it,:) = rad2deg(angle( ...
        complexSpectrum(idx,:,2).*conj(complexSpectrum(idx,:,1))));
end

colors = [0 0 0;0.05 0.35 0.85;0.85 0.15 0.10];
lowMask = frequency<=200;
meshMask = frequency>=120 & frequency<=200;
amplitudeScale = 1e3; % m/s^2 -> mm/s^2 for readable low-frequency axes
fig = figure('Visible','off','Color','w','Position',[80 80 1160 700]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for is = 1:2
    % Complete 0-200 Hz view.
    nexttile((is-1)*2+1);
    curveHandle = gobjects(1,nCase);
    for ic = 1:nCase
        curveHandle(ic) = plot(frequency(lowMask), ...
            amplitudeScale*spectrum(lowMask,ic,is), ...
            'Color',colors(ic,:),'LineWidth',1.0); hold on;
    end
    for harmonic = 1:8
        fFault = harmonic*p.kin.f_sun_fault;
        if fFault<=200 && abs(fFault-p.kin.f_mesh)>1e-9
            xline(fFault,':','Color',[0.55 0.55 0.55], ...
                'HandleVisibility','off');
        end
    end
    xline(p.kin.f_mesh,'--','Color',[0.75 0 0], ...
        'LineWidth',1.2,'HandleVisibility','off');
    xlim([0 200]);
    lowMaximum = amplitudeScale*max(spectrum(lowMask,:,is),[],'all');
    ylim([0 1.08*lowMaximum]);
    xlabel('Frequency (Hz)'); ylabel('Single-sided amplitude (mm/s^2)');
    title(sprintf('Sensor %d deg: 0-200 Hz', ...
        responses{1}.sensorAnglesDeg(is)));
    legend(curveHandle,caseNames,'Location','northwest'); box off;

    % Enlarged mesh band.  Small planet load/path asymmetry restores the
    % fm +/- fc and fm +/- 2*fc orders that cancel in the ideal symmetric
    % case; three-planet passage produces fm +/- 3*fc, and the broken sun
    % tooth produces the wider fm +/- fsf pair.
    nexttile((is-1)*2+2);
    for ic = 1:nCase
        plot(frequency(meshMask),amplitudeScale*spectrum(meshMask,ic,is), ...
            'Color',colors(ic,:),'LineWidth',1.0); hold on;
    end
    xline(p.kin.f_mesh-p.kin.f_sun_fault,':k','144 Hz', ...
        'LabelVerticalAlignment','middle','HandleVisibility','off');
    xline(p.kin.f_mesh-p.model.nPlanet*p.kin.f_c,'-.', ...
        'Color',[0.25 0.45 0.75],'LineWidth',1.0, ...
        'Label','162 Hz path','LabelVerticalAlignment','middle', ...
        'HandleVisibility','off');
    xline(p.kin.f_mesh-p.kin.f_c,'--', ...
        'Color',[0.10 0.55 0.20],'LineWidth',1.0, ...
        'Label','166 Hz: fm-fc','LabelVerticalAlignment','bottom', ...
        'HandleVisibility','off');
    xline(p.kin.f_mesh,'--','168 Hz mesh','Color',[0.75 0 0],'LineWidth',1.2, ...
        'LabelVerticalAlignment','middle', ...
        'HandleVisibility','off');
    xline(p.kin.f_mesh+p.kin.f_c,'--', ...
        'Color',[0.10 0.55 0.20],'LineWidth',1.0, ...
        'Label','170 Hz: fm+fc','LabelVerticalAlignment','top', ...
        'HandleVisibility','off');
    xline(p.kin.f_mesh+p.model.nPlanet*p.kin.f_c,'-.', ...
        'Color',[0.25 0.45 0.75],'LineWidth',1.0, ...
        'Label','174 Hz path','LabelVerticalAlignment','middle', ...
        'HandleVisibility','off');
    xline(p.kin.f_mesh+p.kin.f_sun_fault,':k','192 Hz', ...
        'LabelVerticalAlignment','middle','HandleVisibility','off');
    xlim([120 200]);
    meshMaximum = amplitudeScale*max(spectrum(meshMask,:,is),[],'all');
    ylim([0 1.08*meshMaximum]);
    xlabel('Frequency (Hz)'); ylabel('Single-sided amplitude (mm/s^2)');
    title(sprintf('Sensor %d deg: mesh-frequency detail', ...
        responses{1}.sensorAnglesDeg(is)));
    box off;
end
exportgraphics(fig,fullfile('results','sensor_response_spectrum_0_200Hz.png'), ...
    'Resolution',180);
close(fig);

report.frequencyHz = targetFrequency;
report.names = targetNames;
report.amplitude = targetAmplitude;
report.phaseDifference90Minus0Deg = targetPhaseDifferenceDeg;
report.caseNames = caseNames;
save(fullfile('results','low_frequency_spectrum_report.mat'),'report');

nTarget = numel(targetFrequency);
phaseTable = table('Size',[nTarget*nCase,6], ...
    'VariableTypes',{'string','string','double','double','double','double'}, ...
    'VariableNames',{'case_name','component','frequency_hz', ...
    'amplitude_0deg_ms2','amplitude_90deg_ms2','phase_90_minus_0_deg'});
row = 0;
for ic = 1:nCase
    for it = 1:nTarget
        row = row+1;
        phaseTable.case_name(row) = string(caseNames{ic});
        phaseTable.component(row) = string(targetNames{it});
        phaseTable.frequency_hz(row) = targetFrequency(it);
        phaseTable.amplitude_0deg_ms2(row) = targetAmplitude(it,ic,1);
        phaseTable.amplitude_90deg_ms2(row) = targetAmplitude(it,ic,2);
        phaseTable.phase_90_minus_0_deg(row) = ...
            targetPhaseDifferenceDeg(it,ic);
    end
end
writetable(phaseTable,fullfile('results','lpm_phase_landmarks.csv'));

fprintf('Low-frequency spectrum generated.\n');
for is = 1:2
    fprintf(['Sensor %d deg amplitudes at ', ...
        '144/162/164/166/168/170/172/174/192 Hz:\n'], ...
        responses{1}.sensorAnglesDeg(is));
    for ic = 1:nCase
        fprintf(['  %-13s %.3e / %.3e / %.3e / %.3e / %.3e / ', ...
            '%.3e / %.3e / %.3e / %.3e m/s^2\n'],caseNames{ic}, ...
            targetAmplitude(6,ic,is),targetAmplitude(7,ic,is), ...
            targetAmplitude(8,ic,is),targetAmplitude(9,ic,is), ...
            targetAmplitude(10,ic,is),targetAmplitude(11,ic,is), ...
            targetAmplitude(12,ic,is),targetAmplitude(13,ic,is), ...
            targetAmplitude(14,ic,is));
    end
end
fprintf(['Healthy phase 90-0 at 162/166/168/170/174 Hz: ', ...
    '%.2f / %.2f / %.2f / %.2f / %.2f deg\n'], ...
    targetPhaseDifferenceDeg(7,1),targetPhaseDifferenceDeg(9,1), ...
    targetPhaseDifferenceDeg(10,1),targetPhaseDifferenceDeg(11,1), ...
    targetPhaseDifferenceDeg(13,1));
end
