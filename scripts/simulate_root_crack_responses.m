function output = simulate_root_crack_responses(inputTorqueNm)
%SIMULATE_ROOT_CRACK_RESPONSES One-second loaded sensor root-crack response.

arguments
    inputTorqueNm (1,1) double {mustBeNonnegative} = 20
end

setup_paths;
p = pg_parameters;
p.operating.inputTorque = inputTorqueNm;
if inputTorqueNm > 0
    p.operating.outputTorque = inputTorqueNm*abs(p.kin.omega_s/p.kin.omega_c);
else
    p.operating.outputTorque = 0;
end
% Keep torque ripple separate from the crack mechanism until an encoder- or
% torque-derived amplitude is available.
p.operating.torqueRippleFraction = 0;
% The loaded workflow uses the same provisional 5-um static transmission
% error as loaded-v2. This preserves the physical 168-Hz carrier while the
% crack-induced preload modulation creates smaller surrounding sidebands.
p.mesh.sunPlanet.teAmplitude = 5.0e-6;
p.mesh.ringPlanet.teAmplitude = 5.0e-6;
p.model.name = sprintf('ZLS160 loaded root-crack LPM, Tin=%.4g Nm', ...
    inputTorqueNm);
torqueTag = strrep(sprintf('%.4g',inputTorqueNm),'.','p');
if inputTorqueNm > 0
    fileStem = sprintf('root_crack_loaded%sNm_corrected',torqueTag);
else
    fileStem = 'root_crack_noload_corrected';
end
discardTime = 0.10;
recordDuration = 1.00;
duration = discardTime+recordDuration;
paperDepthsMm = [0,0.5,1.5,2.5];
depthsMm = paperDepthsMm*(1e3*p.gear.module/2.0);
faceCoverage = 0.5;
cases = {pg_fault_case('none',0,p), ...
    pg_root_crack_case(depthsMm(2),p,'FaceCoverage',faceCoverage), ...
    pg_root_crack_case(depthsMm(3),p,'FaceCoverage',faceCoverage), ...
    pg_root_crack_case(depthsMm(4),p,'FaceCoverage',faceCoverage)};
caseNames = {'Healthy',sprintf('Root crack %.3f mm',depthsMm(2)), ...
    sprintf('Root crack %.3f mm',depthsMm(3)), ...
    sprintf('Root crack %.3f mm',depthsMm(4))};

responses = cell(size(cases));
for ic = 1:numel(cases)
    fprintf('Simulating %s ...\n',caseNames{ic});
    response = pg_simulate_lpm_response(duration,p,cases{ic});
    keep = response.time>=discardTime;
    response.time = response.time(keep)-discardTime;
    response.acceleration = response.acceleration(keep,:);
    response.acceleration = response.acceleration-mean(response.acceleration,1);
    assert(all(isfinite(response.acceleration),'all'));
    responses{ic} = response;
end

time = responses{1}.time;
n = numel(time);
nCase = numel(cases);
nSensor = numel(p.sensor.angles);
window = 0.5-0.5*cos(2*pi*(0:n-1)'/(n-1));
frequency = (0:floor(n/2))'*p.model.fs/n;
spectrum = zeros(numel(frequency),nCase,nSensor);
rmsAcceleration = zeros(nCase,nSensor);
crestFactor = zeros(nCase,nSensor);

for ic = 1:nCase
    for is = 1:nSensor
        x = responses{ic}.acceleration(:,is);
        rmsAcceleration(ic,is) = rms(x);
        crestFactor(ic,is) = max(abs(x))/rmsAcceleration(ic,is);
        X = fft(x.*window);
        spectrum(:,ic,is) = 2*abs(X(1:numel(frequency)))/sum(window);
    end
end

% Resonance-band envelope makes the 24-Hz three-planet fault repetition
% directly comparable across crack depths.
frequencyFull = (0:n-1)'*p.model.fs/n;
envelopeBand = [1700,2100];
bandMask = (frequencyFull>=envelopeBand(1) & frequencyFull<=envelopeBand(2)) | ...
    (frequencyFull>=p.model.fs-envelopeBand(2) & ...
     frequencyFull<=p.model.fs-envelopeBand(1));
envelopeSpectrum = zeros(numel(frequency),nCase,nSensor);
for ic = 1:nCase
    for is = 1:nSensor
        bandSignal = real(ifft(fft(responses{ic}.acceleration(:,is)).*bandMask));
        envelope = abs(hilbert(bandSignal));
        envelope = envelope-mean(envelope);
        E = fft(envelope.*window);
        envelopeSpectrum(:,ic,is) = ...
            2*abs(E(1:numel(frequency)))/sum(window);
    end
end

landmarksHz = [p.kin.f_sun_fault,p.kin.f_mesh-p.kin.f_sun_fault, ...
    p.kin.f_mesh,p.kin.f_mesh+p.kin.f_sun_fault];
landmarkAmplitude = zeros(nCase,nSensor,numel(landmarksHz));
faultEnvelopeAmplitude = zeros(nCase,nSensor);
for il = 1:numel(landmarksHz)
    [~,idx] = min(abs(frequency-landmarksHz(il)));
    landmarkAmplitude(:,:,il) = reshape(spectrum(idx,:,:),nCase,nSensor);
end
[~,idxFault] = min(abs(frequency-p.kin.f_sun_fault));
faultEnvelopeAmplitude(:,:) = reshape(envelopeSpectrum(idxFault,:,:), ...
    nCase,nSensor);

% Direct acceptance test requested for the paper figures: the 168-Hz mesh
% carrier must exceed every 6-Hz-grid modulation line in 0-240 Hz.
sidebandGridHz = (6:6:240).';
sidebandGridHz(abs(sidebandGridHz-p.kin.f_mesh)<1e-9) = [];
maximumGridSideband = zeros(nCase,nSensor);
maximumGridSidebandHz = zeros(nCase,nSensor);
meshDominanceMarginDb = zeros(nCase,nSensor);
for ic = 1:nCase
    for is = 1:nSensor
        gridAmplitude = zeros(size(sidebandGridHz));
        for ig = 1:numel(sidebandGridHz)
            [~,idx] = min(abs(frequency-sidebandGridHz(ig)));
            gridAmplitude(ig) = spectrum(idx,ic,is);
        end
        [maximumGridSideband(ic,is),idxMax] = max(gridAmplitude);
        maximumGridSidebandHz(ic,is) = sidebandGridHz(idxMax);
        meshAmplitude = landmarkAmplitude(ic,is,3);
        meshDominanceMarginDb(ic,is) = 20*log10( ...
            meshAmplitude/maximumGridSideband(ic,is));
    end
end
assert(all(meshDominanceMarginDb(:)>0), ...
    'The 168-Hz mesh carrier must exceed every modulation line.');

if ~exist('results','dir'); mkdir('results'); end
metrics = table(repmat(inputTorqueNm,nCase,1), ...
    repmat(p.operating.outputTorque,nCase,1),paperDepthsMm.',depthsMm.', ...
    [0;repmat(faceCoverage,nCase-1,1)], ...
    rmsAcceleration(:,1),rmsAcceleration(:,2), ...
    crestFactor(:,1),crestFactor(:,2),landmarkAmplitude(:,1,2), ...
    landmarkAmplitude(:,1,3),landmarkAmplitude(:,1,4), ...
    faultEnvelopeAmplitude(:,1),faultEnvelopeAmplitude(:,2), ...
    maximumGridSideband(:,1),maximumGridSidebandHz(:,1), ...
    meshDominanceMarginDb(:,1), ...
    'VariableNames',{'input_torque_Nm','carrier_load_Nm', ...
    'paper_crack_depth_mm','scaled_crack_depth_mm','face_coverage', ...
    'rms_0deg_ms2','rms_90deg_ms2', ...
    'crest_0deg','crest_90deg','lower_sideband_0deg_ms2', ...
    'mesh_0deg_ms2','upper_sideband_0deg_ms2', ...
    'fault_envelope_0deg_ms2','fault_envelope_90deg_ms2', ...
    'max_grid_sideband_0deg_ms2','max_grid_sideband_frequency_Hz', ...
    'mesh_dominance_margin_0deg_dB'});
writetable(metrics,fullfile('results',[fileStem '_sensor_metrics.csv']));

signals = table(time);
for ic = 1:nCase
    tag = sprintf('q%04dum',round(1000*depthsMm(ic)));
    signals.([tag '_0deg_ms2']) = responses{ic}.acceleration(:,1);
    signals.([tag '_90deg_ms2']) = responses{ic}.acceleration(:,2);
end
writetable(signals,fullfile('results',[fileStem '_sensor_signals_1s.csv']));

colors = [0.10 0.10 0.10;0.20 0.45 0.75;0.90 0.55 0.10;0.75 0.15 0.18];
fig = figure('Visible','off','Color','w','Position',[80 50 1180 880]);
tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
inTime = time<=0.20;
inBand = frequency<=240;
for ic = 1:nCase
    nexttile;
    plot(time(inTime),responses{ic}.acceleration(inTime,1), ...
        'Color',colors(ic,:),'LineWidth',0.65);
    xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
    title(sprintf('%s: time response',caseNames{ic})); grid on;

    nexttile;
    plot(frequency(inBand),spectrum(inBand,ic,1), ...
        'Color',colors(ic,:),'LineWidth',0.9); hold on;
    xline(p.kin.f_mesh,'--','Color',[0.45 0.45 0.45]);
    xline(p.kin.f_mesh-p.kin.f_sun_fault,':','Color',[0.55 0.55 0.55]);
    xline(p.kin.f_mesh+p.kin.f_sun_fault,':','Color',[0.55 0.55 0.55]);
    xlabel('Frequency (Hz)'); ylabel('Amplitude (m/s^2)');
    title(sprintf('%s: 0-240 Hz',caseNames{ic})); grid on;
end
sgtitle(sprintf('Input torque %.1f N m; carrier resisting torque %.1f N m', ...
    inputTorqueNm,p.operating.outputTorque));
exportgraphics(fig,fullfile('results',[fileStem '_sensor_response.png']), ...
    'Resolution',200);
close(fig);

fig = figure('Visible','off','Color','w','Position',[100 100 1040 620]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
for is = 1:nSensor
    nexttile;
    for ic = 1:nCase
        plot(frequency(frequency<=120), ...
            envelopeSpectrum(frequency<=120,ic,is), ...
            'Color',colors(ic,:),'LineWidth',1.0); hold on;
    end
    for harmonic = 1:5
        xline(harmonic*p.kin.f_sun_fault,':','Color',[0.55 0.55 0.55], ...
            'HandleVisibility','off');
    end
    xlabel('Frequency (Hz)'); ylabel('Envelope amplitude (m/s^2)');
    title(sprintf('Sensor %g deg, %.1f-%.1f kHz envelope', ...
        rad2deg(p.sensor.angles(is)),envelopeBand(1)/1000,envelopeBand(2)/1000));
    legend(caseNames,'Location','best'); grid on;
end
exportgraphics(fig,fullfile('results',[fileStem '_sensor_envelope.png']), ...
    'Resolution',200);
close(fig);

output.caseNames = caseNames;
output.depthsMm = depthsMm;
output.paperDepthsMm = paperDepthsMm;
output.faceCoverage = faceCoverage;
output.inputTorqueNm = inputTorqueNm;
output.carrierLoadNm = p.operating.outputTorque;
output.responses = responses;
output.metrics = metrics;
output.frequency = frequency;
output.spectrum = spectrum;
output.envelopeSpectrum = envelopeSpectrum;
output.landmarksHz = landmarksHz;
output.maximumGridSideband = maximumGridSideband;
output.maximumGridSidebandHz = maximumGridSidebandHz;
output.meshDominanceMarginDb = meshDominanceMarginDb;
save(fullfile('results',[fileStem '_sensor_responses.mat']),'output','-v7.3');

fprintf('Root-crack sensor-response simulation PASSED\n');
fprintf('  input/carrier torque     : %.3f / %.3f N m\n', ...
    inputTorqueNm,p.operating.outputTorque);
fprintf('  TE amplitude             : %.3f um\n', ...
    1e6*p.mesh.sunPlanet.teAmplitude);
fprintf('  scaled q (mm)            : %.3f / %.3f / %.3f\n', ...
    depthsMm(2:4));
fprintf('  mesh dominance 0 deg (dB): %.2f / %.2f / %.2f / %.2f\n', ...
    meshDominanceMarginDb(:,1));
disp(metrics);
end
