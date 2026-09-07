function report = validate_loaded_v2_response
%VALIDATE_LOADED_V2_RESPONSE Mechanism checks for the generated v2 signals.

loaded = load(fullfile('results','loaded_v2_response_cases.mat'),'output');
output = loaded.output;
p = output.parameters;
nCase = numel(output.responses);
nSensor = 2;

minimumCompression = zeros(nCase,1);
meanCompression = zeros(nCase,1);
channelCorrelation = zeros(nCase,1);
contactLossFraction = zeros(nCase,1);
faultOrderAmplitude = zeros(nCase,nSensor,3);
directFaultAmplitude = zeros(nCase,nSensor);
meshCenterAmplitude = zeros(nCase,nSensor);
maximumMeshSidebandAmplitude = zeros(nCase,nSensor);
meshSidebandRatio = zeros(nCase,nSensor);
faultOrders = [1,2,3];

frequency = output.frequency;
[~,faultIndex] = min(abs(frequency-p.kin.f_sun_fault));
[~,meshIndex] = min(abs(frequency-p.kin.f_mesh));
meshBand = frequency>=140 & frequency<=200;
excludeCenter = abs(frequency-p.kin.f_mesh)>1;

for ic = 1:nCase
    compression = -output.responses{ic}.meshForceSunPlanet;
    minimumCompression(ic) = min(compression,[],'all');
    meanCompression(ic) = mean(compression,'all');
    contactLossFraction(ic) = mean(compression<=0,'all');
    correlationMatrix = corrcoef(output.responses{ic}.acceleration(:,1), ...
        output.responses{ic}.acceleration(:,2));
    channelCorrelation(ic) = correlationMatrix(1,2);

    x = output.responses{ic}.acceleration;
    fs = output.responses{ic}.fs;
    n = size(x,1);
    frequencyFull = (0:n-1)'*fs/n;
    band = [1700 2050];
    bandMask = (frequencyFull>=band(1) & frequencyFull<=band(2)) | ...
        (frequencyFull>=fs-band(2) & frequencyFull<=fs-band(1));
    window = 0.5-0.5*cos(2*pi*(0:n-1)'/(n-1));
    for is = 1:nSensor
        rawAmplitude = output.spectrum(:,ic,is);
        directFaultAmplitude(ic,is) = rawAmplitude(faultIndex);
        meshCenterAmplitude(ic,is) = rawAmplitude(meshIndex);
        maximumMeshSidebandAmplitude(ic,is) = ...
            max(rawAmplitude(meshBand & excludeCenter));
        meshSidebandRatio(ic,is) = maximumMeshSidebandAmplitude(ic,is)/ ...
            max(meshCenterAmplitude(ic,is),eps);

        bandSignal = real(ifft(fft(x(:,is)).*bandMask));
        envelope = abs(hilbert(bandSignal));
        envelope = envelope-mean(envelope);
        E = fft(envelope.*window);
        envelopeAmplitude = 2*abs(E(1:numel(frequency)))/sum(window);
        for io = 1:numel(faultOrders)
            target = faultOrders(io)*p.kin.f_sun_fault;
            [~,index] = min(abs(frequency-target));
            faultOrderAmplitude(ic,is,io) = envelopeAmplitude(index);
        end
    end
end

assert(all(contactLossFraction(1:2)==0), ...
    'Healthy or 25-percent case lost compressive contact.');
assert(contactLossFraction(3)<1e-3, ...
    'Severe-case contact loss is too large for the linear-contact model.');
assert(all(cellfun(@(r) abs(r.fs-51200)<eps(51200),output.responses)), ...
    'A response was not generated at 51.2 kHz.');
assert(all(diff(output.rmsAcceleration(:,1))>0), ...
    'Sensor-0 RMS does not increase monotonically with severity.');
assert(all(diff(faultOrderAmplitude(:,1,1))>0), ...
    'The 24-Hz envelope indicator does not increase with severity.');
assert(all(meshSidebandRatio(:,1)<=1.6), ...
    'A mesh-band sideband exceeds the 168-Hz center by more than 60 percent.');

report.minimumCompressionN = minimumCompression;
report.meanCompressionN = meanCompression;
report.channelCorrelation = channelCorrelation;
report.contactLossFraction = contactLossFraction;
report.faultOrders = faultOrders;
report.faultOrderAmplitude = faultOrderAmplitude;
report.directFaultAmplitude = directFaultAmplitude;
report.meshCenterAmplitude = meshCenterAmplitude;
report.maximumMeshSidebandAmplitude = maximumMeshSidebandAmplitude;
report.meshSidebandRatio = meshSidebandRatio;

summary = table(string(output.caseNames(:)),minimumCompression, ...
    meanCompression,contactLossFraction,channelCorrelation, ...
    squeeze(faultOrderAmplitude(:,1,1)), ...
    squeeze(faultOrderAmplitude(:,1,2)), ...
    squeeze(faultOrderAmplitude(:,1,3)), ...
    directFaultAmplitude(:,1),meshCenterAmplitude(:,1), ...
    maximumMeshSidebandAmplitude(:,1),meshSidebandRatio(:,1), ...
    'VariableNames',{'case_name','minimum_compression_N', ...
    'mean_compression_N','contact_loss_fraction', ...
    'zero_lag_channel_correlation', ...
    'envelope_24Hz_sensor0_ms2','envelope_48Hz_sensor0_ms2', ...
    'envelope_72Hz_sensor0_ms2','raw_24Hz_sensor0_ms2', ...
    'mesh_168Hz_sensor0_ms2','max_mesh_sideband_sensor0_ms2', ...
    'max_sideband_to_mesh_ratio'});
writetable(summary,fullfile('results','loaded_v2_validation.csv'));

fprintf('Loaded-v2 validation PASSED.\n');
for ic = 1:nCase
    fprintf(['  %-24s min/mean compression %.1f/%.1f N; loss %.5f%%; ' ...
        'channel corr %.3f; raw 24/mesh/max-side %.4g/%.4g/%.4g; ' ...
        'side/mesh %.3f; envelope 24/48/72 Hz %.3f/%.3f/%.3f m/s^2\n'], ...
        output.caseNames{ic},minimumCompression(ic),meanCompression(ic), ...
        100*contactLossFraction(ic),channelCorrelation(ic), ...
        directFaultAmplitude(ic,1),meshCenterAmplitude(ic,1), ...
        maximumMeshSidebandAmplitude(ic,1),meshSidebandRatio(ic,1), ...
        faultOrderAmplitude(ic,1,1), ...
        faultOrderAmplitude(ic,1,2),faultOrderAmplitude(ic,1,3));
end
end
