function report = analyze_loaded_path_am_pm
%ANALYZE_LOADED_PATH_AM_PM Ring-path AM/PM study for the loaded v2.1 model.

setup_paths;
p = pg_parameters_loaded_v2;
model = pg_ring_modal_model(p);

% The 11th mesh harmonic lies in the first measured elastic-ring band.
analysisFrequency = [p.kin.f_mesh,11*p.kin.f_mesh];
analysisName = {'168-Hz mesh center','1848-Hz resonance carrier'};
phi = linspace(0,2*pi,721);
pointField = pg_path_field(analysisFrequency,phi,p,model,ones(1,3));
singlePlanetField = pg_path_field(analysisFrequency,phi,p,model,[1,0,0]);

nFrequency = numel(analysisFrequency);
nSensor = 2;
amplitudeCv = zeros(nFrequency,nSensor);
amplitudeDepth = zeros(nFrequency,nSensor);
phaseSpanDeg = zeros(nFrequency,nSensor);
dualPhaseSpanDeg = zeros(nFrequency,1);

for jf = 1:nFrequency
    for is = 1:nSensor
        h = squeeze(pointField.HCombined(jf,:,is));
        a = abs(h);
        amplitudeCv(jf,is) = std(a)/mean(a);
        amplitudeDepth(jf,is) = (max(a)-min(a))/(max(a)+min(a)+eps);
        valid = a>0.05*max(a);
        phase = unwrap(angle(h));
        phaseSpanDeg(jf,is) = rad2deg(range(phase(valid)));
    end
    phaseDifference = pointField.phaseDifference(jf,:);
    valid = pointField.validMask(jf,:);
    dualPhaseSpanDeg(jf) = rad2deg(circularSpan(phaseDifference(valid)));
end

frequencyColumn = repelem(analysisFrequency(:),nSensor);
sensorColumn = repmat([0;90],nFrequency,1);
metrics = table(frequencyColumn,sensorColumn,reshape(amplitudeCv.',[],1), ...
    reshape(amplitudeDepth.',[],1),reshape(phaseSpanDeg.',[],1), ...
    repelem(dualPhaseSpanDeg,nSensor), ...
    'VariableNames',{'frequency_hz','sensor_angle_deg','amplitude_cv', ...
    'amplitude_depth','path_phase_span_deg','dual_phase_span_deg'});
writetable(metrics,fullfile('results','loaded_path_am_pm_metrics.csv'));

% Angle-domain path amplitude and phase modulation.
fig = figure('Visible','off','Color','w','Position',[70 50 1180 820]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for jf = 1:nFrequency
    nexttile;
    magnitude = squeeze(abs(pointField.HCombined(jf,:,:)));
    magnitude = magnitude/max(magnitude,[],'all');
    plot(rad2deg(phi),magnitude(:,1),'b','LineWidth',1.1); hold on;
    plot(rad2deg(phi),magnitude(:,2),'r','LineWidth',1.1);
    xlim([0 360]); ylim([0 1.05]); grid on; box off;
    xlabel('Carrier angle (deg)'); ylabel('Normalized path magnitude');
    title(sprintf('%s: amplitude modulation',analysisName{jf}));
    legend('0 deg sensor','90 deg sensor','Location','best');

    nexttile;
    phaseDifference = rad2deg(pointField.phaseDifference(jf,:));
    valid = pointField.validMask(jf,:);
    plot(rad2deg(phi(valid)),phaseDifference(valid),'k.','MarkerSize',5);
    xlim([0 360]); ylim([-180 180]); grid on; box off;
    xlabel('Carrier angle (deg)'); ylabel('H90-H0 phase (deg)');
    title(sprintf('%s: phase modulation',analysisName{jf}));
end
exportgraphics(fig,fullfile('results','loaded_path_am_pm.png'), ...
    'Resolution',190);
close(fig);

% A fault pulse is localized at one mesh, so retain the single-planet path
% rather than assigning the coherent three-planet phase to every event.
singleDualPhaseSpanDeg = zeros(nFrequency,1);
fig = figure('Visible','off','Color','w','Position',[70 50 1180 820]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for jf = 1:nFrequency
    nexttile;
    magnitude = squeeze(abs(singlePlanetField.HCombined(jf,:,:)));
    magnitude = magnitude/max(magnitude,[],'all');
    plot(rad2deg(phi),magnitude(:,1),'b','LineWidth',1.1); hold on;
    plot(rad2deg(phi),magnitude(:,2),'r','LineWidth',1.1);
    xlim([0 360]); ylim([0 1.05]); grid on; box off;
    xlabel('Carrier angle (deg)'); ylabel('Normalized path magnitude');
    title(sprintf('%s: planet 1 path AM',analysisName{jf}));
    legend('0 deg sensor','90 deg sensor','Location','best');

    nexttile;
    h0 = squeeze(singlePlanetField.HCombined(jf,:,1));
    h90 = squeeze(singlePlanetField.HCombined(jf,:,2));
    [singleReliability,singleValid] = pathReliability(h0,h90);
    phaseDifference = angle(h90.*conj(h0));
    singleDualPhaseSpanDeg(jf) = ...
        rad2deg(circularSpan(phaseDifference(singleValid)));
    scatter(rad2deg(phi(singleValid)),rad2deg(phaseDifference(singleValid)), ...
        8,singleReliability(singleValid),'filled');
    xlim([0 360]); ylim([-180 180]); grid on; box off;
    xlabel('Carrier angle (deg)'); ylabel('H90-H0 phase (deg)');
    title(sprintf('%s: planet 1 path PM',analysisName{jf}));
    colorbar;
end
exportgraphics(fig,fullfile('results', ...
    'loaded_single_planet_path_am_pm.png'),'Resolution',190);
close(fig);

% Fourier-series decomposition of a periodic complex path. Its coefficients
% are path-induced carrier-order sidebands f0+k*fc.
orders = -12:12;
loaded = load(fullfile('results','loaded_v2_response_cases.mat'),'output');
output = loaded.output;
fullCoefficient = zeros(nFrequency,numel(orders));
amCoefficient = zeros(size(fullCoefficient));
pmCoefficient = zeros(size(fullCoefficient));
measuredAmplitude = zeros(size(fullCoefficient));
sidebandCorrelation = zeros(nFrequency,1);
sidebandNrmse = zeros(nFrequency,1);

fig = figure('Visible','off','Color','w','Position',[80 80 1120 700]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
for jf = 1:nFrequency
    h = squeeze(pointField.HCombined(jf,1:end-1,1));
    referencePhase = angle(mean(h));
    hAmplitudeOnly = abs(h).*exp(1i*referencePhase);
    hPhaseOnly = mean(abs(h))*exp(1i*angle(h));
    fullCoefficient(jf,:) = fourierOrders(h,orders);
    amCoefficient(jf,:) = fourierOrders(hAmplitudeOnly,orders);
    pmCoefficient(jf,:) = fourierOrders(hPhaseOnly,orders);

    for io = 1:numel(orders)
        target = analysisFrequency(jf)+orders(io)*p.kin.f_c;
        [~,index] = min(abs(output.frequency-target));
        measuredAmplitude(jf,io) = output.spectrum(index,1,1);
    end

    scalePath = max(abs(fullCoefficient(jf,:)));
    scaleMeasured = max(measuredAmplitude(jf,:));
    nexttile;
    bar(orders-0.24,abs(fullCoefficient(jf,:))/scalePath,0.23, ...
        'FaceColor',[0.15 0.15 0.15],'EdgeColor','none'); hold on;
    bar(orders,abs(amCoefficient(jf,:))/scalePath,0.23, ...
        'FaceColor',[0.10 0.45 0.85],'EdgeColor','none');
    bar(orders+0.24,abs(pmCoefficient(jf,:))/scalePath,0.23, ...
        'FaceColor',[0.85 0.25 0.10],'EdgeColor','none');
    plot(orders,measuredAmplitude(jf,:)/max(scaleMeasured,eps), ...
        'ko-','LineWidth',1.0,'MarkerSize',4,'MarkerFaceColor','w');
    xlim([orders(1)-0.7 orders(end)+0.7]); ylim([0 1.08]);
    grid on; box off;
    xlabel('Carrier sideband order k'); ylabel('Normalized amplitude');
    title(sprintf('%s: f_0+k f_c sidebands',analysisName{jf}));
    legend('Full complex path','AM only','PM only', ...
        'Healthy simulation','Location','best');

    normalizedPath = abs(fullCoefficient(jf,:))/scalePath;
    normalizedMeasured = measuredAmplitude(jf,:)/max(scaleMeasured,eps);
    correlationMatrix = corrcoef(normalizedPath,normalizedMeasured);
    sidebandCorrelation(jf) = correlationMatrix(1,2);
    sidebandNrmse(jf) = sqrt(mean((normalizedPath-normalizedMeasured).^2));
end
exportgraphics(fig,fullfile('results', ...
    'loaded_path_sideband_decomposition.png'),'Resolution',190);
close(fig);

sidebandTable = table;
for jf = 1:nFrequency
    block = table(repmat(analysisFrequency(jf),numel(orders),1), ...
        orders(:),abs(fullCoefficient(jf,:)).', ...
        abs(amCoefficient(jf,:)).',abs(pmCoefficient(jf,:)).', ...
        measuredAmplitude(jf,:).', ...
        'VariableNames',{'carrier_frequency_hz','carrier_order', ...
        'full_path_coefficient','am_only_coefficient', ...
        'pm_only_coefficient','healthy_simulation_amplitude'});
    sidebandTable = [sidebandTable;block]; %#ok<AGROW>
end
writetable(sidebandTable,fullfile('results', ...
    'loaded_path_sideband_decomposition.csv'));

% Frequency-angle maps for the mesh band and first resonance band.
lowFrequency = (140:1:200).';
resonanceFrequency = (1700:5:2050).';
mapFrequency = [lowFrequency;resonanceFrequency];
phiMap = linspace(0,2*pi,181);
mapField = pg_path_field(mapFrequency,phiMap,p,model,ones(1,3));
fig = figure('Visible','off','Color','w','Position',[50 40 1180 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
plotBandMap(lowFrequency,1:numel(lowFrequency),'140-200 Hz mesh band');
offset = numel(lowFrequency);
plotBandMap(resonanceFrequency,offset+(1:numel(resonanceFrequency)), ...
    '1.70-2.05 kHz resonance band');
exportgraphics(fig,fullfile('results','loaded_path_truth_maps.png'), ...
    'Resolution',190);
close(fig);

% Encoder-synchronous phase truth for every possible dominant planet path.
encoderAngle = 2*pi*(0:p.operating.encoderPpr-1)/p.operating.encoderPpr;
lookupField = pg_path_field(analysisFrequency,encoderAngle,p,model,ones(1,3));
lookupTable = table;
for jf = 1:nFrequency
    for ip = 1:p.model.nPlanet
        h0 = squeeze(lookupField.H(jf,:,1,ip));
        h90 = squeeze(lookupField.H(jf,:,2,ip));
        [reliability,valid] = pathReliability(h0,h90);
        phaseDifference = angle(h90.*conj(h0));
        block = table(repmat(analysisFrequency(jf),numel(encoderAngle),1), ...
            repmat(ip,numel(encoderAngle),1),rad2deg(encoderAngle(:)), ...
            abs(h0(:)),abs(h90(:)),rad2deg(angle(h0(:))), ...
            rad2deg(angle(h90(:))),rad2deg(phaseDifference(:)), ...
            reliability(:),valid(:), ...
            'VariableNames',{'frequency_hz','planet_index', ...
            'carrier_angle_deg','H0_magnitude','H90_magnitude', ...
            'H0_phase_deg','H90_phase_deg','phase_difference_deg', ...
            'reliability','valid'});
        lookupTable = [lookupTable;block]; %#ok<AGROW>
    end
end
writetable(lookupTable,fullfile('results','loaded_path_phase_lookup.csv'));
lookupValidFraction = zeros(nFrequency,1);
lookupMedianReliability = zeros(nFrequency,1);
for jf = 1:nFrequency
    useFrequency = lookupTable.frequency_hz==analysisFrequency(jf);
    lookupValidFraction(jf) = mean(lookupTable.valid(useFrequency));
    lookupMedianReliability(jf) = ...
        median(lookupTable.reliability(useFrequency));
end

validationTable = table(analysisFrequency(:),sidebandCorrelation, ...
    sidebandNrmse,lookupValidFraction,lookupMedianReliability, ...
    dualPhaseSpanDeg,singleDualPhaseSpanDeg, ...
    'VariableNames',{'frequency_hz','path_simulation_sideband_correlation', ...
    'normalized_sideband_rmse','lookup_valid_fraction', ...
    'lookup_median_reliability','combined_dual_phase_span_deg', ...
    'single_planet_dual_phase_span_deg'});
writetable(validationTable,fullfile('results', ...
    'loaded_path_am_pm_validation.csv'));

report.parameters = p;
report.model = model;
report.pointField = pointField;
report.singlePlanetField = singlePlanetField;
report.mapField = mapField;
report.lookupTable = lookupTable;
report.validationTable = validationTable;
report.metrics = metrics;
report.sidebandTable = sidebandTable;
report.analysisFrequencyHz = analysisFrequency;
report.analysisDefinition = [ ...
    'AM and PM are generated by the carrier-angle-periodic complex ' ...
    'ring path; Fourier order k maps to f0+k*fc.'];
save(fullfile('results','loaded_path_am_pm.mat'),'report','-v7.3');

fprintf('Loaded path AM/PM analysis complete.\n');
for jf = 1:nFrequency
    fprintf(['  %.0f Hz: amplitude CV 0/90 %.3f/%.3f; ' ...
        'dual-phase span %.1f deg\n'],analysisFrequency(jf), ...
        amplitudeCv(jf,1),amplitudeCv(jf,2),dualPhaseSpanDeg(jf));
    fprintf('           planet-1 reliable dual-phase span %.1f deg\n', ...
        singleDualPhaseSpanDeg(jf));
    fprintf(['           sideband corr/NRMSE %.5f/%.4f; ' ...
        'lookup valid %.3f\n'],sidebandCorrelation(jf), ...
        sidebandNrmse(jf),lookupValidFraction(jf));
end

    function plotBandMap(frequency,index,label)
        combined = mapField.HCombined(index,:,:);
        magnitude = sqrt(abs(combined(:,:,1)).*abs(combined(:,:,2)));
        magnitudeDb = 20*log10(magnitude/max(magnitude,[],'all')+eps);
        phase = rad2deg(mapField.phaseDifference(index,:));
        reliable = mapField.validMask(index,:) & magnitudeDb>-45;
        phase(~reliable) = NaN;

        nexttile;
        imagesc(rad2deg(phiMap),frequency,magnitudeDb);
        axis xy; caxis([-45 0]); colorbar;
        xlabel('Carrier angle (deg)'); ylabel('Frequency (Hz)');
        title([label,': dual-path magnitude (dB)']);

        nexttile;
        imagesc(rad2deg(phiMap),frequency,phase);
        axis xy; caxis([-180 180]); colorbar;
        xlabel('Carrier angle (deg)'); ylabel('Frequency (Hz)');
        title([label,': H90-H0 phase (deg)']);
    end
end

function span = circularSpan(angleValue)
angleValue = mod(angleValue(:),2*pi);
angleValue = sort(angleValue(isfinite(angleValue)));
if numel(angleValue) < 2
    span = 0;
    return;
end
gaps = diff([angleValue;angleValue(1)+2*pi]);
span = 2*pi-max(gaps);
end

function [reliability,valid] = pathReliability(h0,h90)
amplitude0 = abs(h0);
amplitude90 = abs(h90);
relative0 = amplitude0/max(amplitude0);
relative90 = amplitude90/max(amplitude90);
balance = 2*amplitude0.*amplitude90./ ...
    (amplitude0.^2+amplitude90.^2+eps);
reliability = balance.*sqrt(relative0.*relative90);
valid = reliability>0.05 & relative0>1e-4 & relative90>1e-4;
end

function coefficient = fourierOrders(periodicPath,orders)
periodicPath = periodicPath(:).';
n = numel(periodicPath);
transform = fft(periodicPath)/n;
coefficient = complex(zeros(size(orders)));
for k = 1:numel(orders)
    index = mod(orders(k),n)+1;
    coefficient(k) = transform(index);
end
end
