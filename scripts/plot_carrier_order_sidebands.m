function report = plot_carrier_order_sidebands
%PLOT_CARRIER_ORDER_SIDEBANDS Show fm+n*fc components on a relative dB scale.

setup_paths;
p = pg_parameters;
loaded = load(fullfile('results','sensor_response_sun_cases_1s.mat'),'output');
output = loaded.output;
responses = output.responses;
caseNames = output.caseNames;
nCase = numel(responses);
n = numel(responses{1}.time);
fs = responses{1}.fs;
window = 0.5-0.5*cos(2*pi*(0:n-1)'/(n-1));
frequency = (0:floor(n/2))'*fs/n;

carrierOrder = -6:6;
targetFrequency = p.kin.f_mesh+carrierOrder*p.kin.f_c;
amplitude = zeros(numel(carrierOrder),nCase,2);
phaseDifferenceDeg = zeros(numel(carrierOrder),nCase);

for ic = 1:nCase
    complexSpectrum = complex(zeros(numel(frequency),2));
    for is = 1:2
        x = responses{ic}.acceleration(:,is);
        x = x-mean(x);
        X = fft(x.*window);
        complexSpectrum(:,is) = 2*X(1:numel(frequency))/sum(window);
    end
    for io = 1:numel(carrierOrder)
        [~,idx] = min(abs(frequency-targetFrequency(io)));
        amplitude(io,ic,:) = abs(complexSpectrum(idx,:));
        phaseDifferenceDeg(io,ic) = rad2deg(angle( ...
            complexSpectrum(idx,2)*conj(complexSpectrum(idx,1))));
    end
end

centerIndex = find(carrierOrder==0,1);
relativeDb = zeros(size(amplitude));
for ic = 1:nCase
    for is = 1:2
        reference = amplitude(centerIndex,ic,is);
        relativeDb(:,ic,is) = 20*log10(max(amplitude(:,ic,is),eps)/reference);
    end
end

colors = [0 0 0;0.05 0.35 0.85;0.85 0.15 0.10];
fig = figure('Visible','off','Color','w','Position',[100 100 1040 620]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
for is = 1:2
    nexttile;
    curveHandle = gobjects(1,nCase);
    for ic = 1:nCase
        curveHandle(ic) = plot(targetFrequency,relativeDb(:,ic,is),'-o', ...
            'Color',colors(ic,:),'LineWidth',1.0,'MarkerSize',4); hold on;
    end
    xline(p.kin.f_mesh,'--','168 Hz mesh','Color',[0.75 0 0], ...
        'LineWidth',1.2,'HandleVisibility','off');
    for io = 1:numel(carrierOrder)
        xline(targetFrequency(io),':','Color',[0.75 0.75 0.75], ...
            'HandleVisibility','off');
    end
    xlim([targetFrequency(1),targetFrequency(end)]);
    ylim([-100 5]);
    xticks(targetFrequency);
    xlabel('Frequency f_m+n f_c (Hz)');
    ylabel('Amplitude relative to 168 Hz (dB)');
    title(sprintf('Sensor %d deg: carrier-order sidebands, f_c = %.0f Hz', ...
        responses{1}.sensorAnglesDeg(is),p.kin.f_c));
    legend(curveHandle,caseNames,'Location','southoutside', ...
        'Orientation','horizontal');
    grid on;
end
exportgraphics(fig,fullfile('results','carrier_order_sidebands.png'), ...
    'Resolution',180);
close(fig);

report.carrierOrder = carrierOrder;
report.frequencyHz = targetFrequency;
report.amplitude = amplitude;
report.relativeDb = relativeDb;
report.phaseDifference90Minus0Deg = phaseDifferenceDeg;
report.caseNames = caseNames;
save(fullfile('results','carrier_order_sidebands.mat'),'report');

nRow = numel(carrierOrder)*nCase;
resultTable = table('Size',[nRow,8], ...
    'VariableTypes',{'string','double','double','double','double', ...
    'double','double','double'}, ...
    'VariableNames',{'case_name','carrier_order_n','frequency_hz', ...
    'amplitude_0deg_ms2','amplitude_90deg_ms2','relative_0deg_db', ...
    'relative_90deg_db','phase_90_minus_0_deg'});
row = 0;
for ic = 1:nCase
    for io = 1:numel(carrierOrder)
        row = row+1;
        resultTable.case_name(row) = string(caseNames{ic});
        resultTable.carrier_order_n(row) = carrierOrder(io);
        resultTable.frequency_hz(row) = targetFrequency(io);
        resultTable.amplitude_0deg_ms2(row) = amplitude(io,ic,1);
        resultTable.amplitude_90deg_ms2(row) = amplitude(io,ic,2);
        resultTable.relative_0deg_db(row) = relativeDb(io,ic,1);
        resultTable.relative_90deg_db(row) = relativeDb(io,ic,2);
        resultTable.phase_90_minus_0_deg(row) = phaseDifferenceDeg(io,ic);
    end
end
writetable(resultTable,fullfile('results','carrier_order_sidebands.csv'));

fprintf('Carrier-order sideband plot generated.\n');
fprintf('Healthy sensor-0 relative levels n=-3..3 (dB):');
selected = carrierOrder>=-3 & carrierOrder<=3;
fprintf(' %.2f',relativeDb(selected,1,1));
fprintf('\n');
end
