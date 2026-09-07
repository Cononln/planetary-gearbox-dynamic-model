function report = plot_mesh_band_carrier_envelope
%PLOT_MESH_BAND_CARRIER_ENVELOPE Demodulate fc orders around the 168-Hz mesh.

setup_paths;
p = pg_parameters;
loaded = load(fullfile('results','sensor_response_sun_cases_1s.mat'),'output');
output = loaded.output;
responses = output.responses;
caseNames = output.caseNames;
nCase = numel(responses);
n = numel(responses{1}.time);
fs = responses{1}.fs;
frequencyFull = (0:n-1)'*fs/n;
frequency = (0:floor(n/2))'*fs/n;
window = 0.5-0.5*cos(2*pi*(0:n-1)'/(n-1));

meshBand = [150,186];
bandMask = (frequencyFull>=meshBand(1) & frequencyFull<=meshBand(2)) | ...
    (frequencyFull>=fs-meshBand(2) & frequencyFull<=fs-meshBand(1));
envelopeSpectrum = zeros(numel(frequency),nCase,2);

for ic = 1:nCase
    for is = 1:2
        x = responses{ic}.acceleration(:,is);
        x = x-mean(x);
        meshBandSignal = real(ifft(fft(x).*bandMask));
        envelope = abs(hilbert(meshBandSignal));
        envelope = envelope-mean(envelope);
        E = fft(envelope.*window);
        envelopeSpectrum(:,ic,is) = ...
            2*abs(E(1:numel(frequency)))/sum(window);
    end
end

orders = 1:6;
orderFrequency = orders*p.kin.f_c;
orderAmplitude = zeros(numel(orders),nCase,2);
for io = 1:numel(orders)
    [~,idx] = min(abs(frequency-orderFrequency(io)));
    orderAmplitude(io,:,:) = envelopeSpectrum(idx,:,:);
end

displayScale = 1e3;
plotMask = frequency<=15;
colors = [0 0 0;0.05 0.35 0.85;0.85 0.15 0.10];
fig = figure('Visible','off','Color','w','Position',[100 100 1040 620]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
for is = 1:2
    nexttile;
    curveHandle = gobjects(1,nCase);
    for ic = 1:nCase
        curveHandle(ic) = plot(frequency(plotMask), ...
            displayScale*envelopeSpectrum(plotMask,ic,is), ...
            'Color',colors(ic,:),'LineWidth',1.1); hold on;
    end
    xline(p.kin.f_c,'--','2 Hz = f_c','Color',[0.10 0.55 0.20], ...
        'LineWidth',1.1,'HandleVisibility','off');
    xline(2*p.kin.f_c,'--','4 Hz = 2f_c','Color',[0.35 0.55 0.20], ...
        'LineWidth',1.0,'HandleVisibility','off');
    xline(3*p.kin.f_c,'--','6 Hz = 3f_c','Color',[0.25 0.45 0.75], ...
        'LineWidth',1.1,'HandleVisibility','off');
    xlim([0 15]);
    xlabel('Envelope frequency (Hz)');
    ylabel('Envelope amplitude (mm/s^2)');
    title(sprintf(['Sensor %d deg: %.0f-%.0f Hz mesh-band envelope, ', ...
        'f_c = %.0f Hz'],responses{1}.sensorAnglesDeg(is), ...
        meshBand(1),meshBand(2),p.kin.f_c));
    legend(curveHandle,caseNames,'Location','northeast');
    grid on;
end
exportgraphics(fig,fullfile('results','mesh_band_carrier_envelope.png'), ...
    'Resolution',180);
close(fig);

report.meshBandHz = meshBand;
report.order = orders;
report.frequencyHz = orderFrequency;
report.amplitude = orderAmplitude;
report.caseNames = caseNames;
save(fullfile('results','mesh_band_carrier_envelope.mat'),'report');

nRow = numel(orders)*nCase;
resultTable = table('Size',[nRow,5], ...
    'VariableTypes',{'string','double','double','double','double'}, ...
    'VariableNames',{'case_name','carrier_order_n','frequency_hz', ...
    'amplitude_0deg_ms2','amplitude_90deg_ms2'});
row = 0;
for ic = 1:nCase
    for io = 1:numel(orders)
        row = row+1;
        resultTable.case_name(row) = string(caseNames{ic});
        resultTable.carrier_order_n(row) = orders(io);
        resultTable.frequency_hz(row) = orderFrequency(io);
        resultTable.amplitude_0deg_ms2(row) = orderAmplitude(io,ic,1);
        resultTable.amplitude_90deg_ms2(row) = orderAmplitude(io,ic,2);
    end
end
writetable(resultTable, ...
    fullfile('results','mesh_band_carrier_envelope.csv'));

fprintf('Mesh-band carrier envelope generated.\n');
fprintf('Healthy sensor-0 fc/2fc/3fc amplitudes: %.3e / %.3e / %.3e m/s^2\n', ...
    orderAmplitude(1,1,1),orderAmplitude(2,1,1),orderAmplitude(3,1,1));
end
