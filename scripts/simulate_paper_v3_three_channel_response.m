function output = simulate_paper_v3_three_channel_response(options)
%SIMULATE_PAPER_V3_THREE_CHANNEL_RESPONSE Generate the paper-facing model.
% The workflow uses the actual 0/120/240-degree sensor geometry and writes
% source fault-event truth alongside the three observed acceleration traces.

arguments
    options.fs (1,1) double {mustBePositive} = 51200
    options.recordDuration (1,1) double {mustBePositive} = 3.0
    options.discardDuration (1,1) double {mustBeNonnegative} = 0.2
    options.includeSpeedVariation (1,1) logical = false
    options.speedVariationFraction (1,1) double {mustBeNonnegative} = 0.015
    options.speedVariationFrequencyHz (1,1) double {mustBeNonnegative} = 0.65
end

setup_paths;
if options.speedVariationFraction>=1
    error('simulate_paper_v3_three_channel_response:SpeedVariation', ...
        'speedVariationFraction must be below one.');
end
p = pg_parameters_paper_v3;
if options.includeSpeedVariation
    p.operating.speedProfile.enabled = true;
    p.operating.speedProfile.fraction = options.speedVariationFraction;
    p.operating.speedProfile.frequencyHz = options.speedVariationFrequencyHz;
end
model = pg_ring_modal_model(p);
faultCases = {pg_fault_case('none',0,p), ...
              pg_fault_case('sun',0.25,p), ...
              pg_fault_case('sun',0.50,p)};
caseNames = {'Healthy','Sun broken tooth 25%','Sun broken tooth 50%'};
caseTokens = {'healthy','sun25','sun50'};
duration = options.recordDuration+options.discardDuration;
nCase = numel(faultCases);
nSensor = numel(p.sensor.locationAngles);
responses = cell(nCase,1);

for ic = 1:nCase
    fprintf('Paper-v3 simulation: %s\n',caseNames{ic});
    response = pg_simulate_loaded_v2_response(duration,p,faultCases{ic}, ...
        model,'fs',options.fs);
    keep = response.time>=options.discardDuration;
    response.time = response.time(keep)-options.discardDuration;
    response.acceleration = response.acceleration(keep,:);
    response.acceleration = response.acceleration-mean(response.acceleration,1);
    response.meshForceSunPlanet = response.meshForceSunPlanet(keep,:);
    response.meshForceRingPlanet = response.meshForceRingPlanet(keep,:);
    response.inputTorque = response.inputTorque(keep,:);
    response.carrierAngleRad = response.carrierAngleRad(keep,:);
    response.carrierSpeedRadS = response.carrierSpeedRadS(keep,:);
    response.meshPhaseRad = response.meshPhaseRad(keep,:);
    response.meshFrequencyHz = response.meshFrequencyHz(keep,:);
    response.faultStiffnessLossSunPlanet = ...
        response.faultStiffnessLossSunPlanet(keep,:);
    response.faultStiffnessLossRingPlanet = ...
        response.faultStiffnessLossRingPlanet(keep,:);
    response.faultEventTruth = pg_fault_event_truth(response.time, ...
        response.faultStiffnessLossSunPlanet+ ...
        response.faultStiffnessLossRingPlanet,p,faultCases{ic});
    assert(all(isfinite(response.acceleration),'all'));
    responses{ic} = response;
end

time = responses{1}.time;
nSample = numel(time);
window = 0.5-0.5*cos(2*pi*(0:nSample-1)'/(nSample-1));
frequency = (0:floor(nSample/2))'*options.fs/nSample;
spectrum = zeros(numel(frequency),nCase,nSensor);
rmsAcceleration = zeros(nCase,nSensor);
crestFactor = zeros(nCase,nSensor);
kurtosisValue = zeros(nCase,nSensor);
for ic = 1:nCase
    x = responses{ic}.acceleration;
    rmsAcceleration(ic,:) = sqrt(mean(x.^2,1));
    crestFactor(ic,:) = max(abs(x),[],1)./max(rmsAcceleration(ic,:),eps);
    variance = mean(x.^2,1);
    kurtosisValue(ic,:) = mean(x.^4,1)./max(variance.^2,eps);
    for is = 1:nSensor
        X = fft(x(:,is).*window);
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
output.speedVariationEnabled = options.includeSpeedVariation;

if options.includeSpeedVariation
    suffix = '_speedvar';
else
    suffix = '_nominal';
end
save(fullfile('results',['paper_v3_three_channel_response',suffix,'.mat']), ...
    'output','-v7.3');

signalMatrix = time;
signalNames = {'time_s'};
for ic = 1:nCase
    for is = 1:nSensor
        signalMatrix(:,end+1) = responses{ic}.acceleration(:,is); %#ok<AGROW>
        signalNames{end+1} = sprintf('%s_sensor_%03ddeg_ms2', ...
            caseTokens{ic},round(responses{ic}.sensorAnglesDeg(is))); %#ok<AGROW>
    end
end
writetable(array2table(signalMatrix,'VariableNames',signalNames), ...
    fullfile('results',['paper_v3_three_channel_signals',suffix,'.csv']));

metrics = table(string(caseNames(:)),'VariableNames',{'case_name'});
for is = 1:nSensor
    angleToken = sprintf('s%03ddeg',round(responses{1}.sensorAnglesDeg(is)));
    metrics.(['rms_',angleToken,'_ms2']) = rmsAcceleration(:,is);
    metrics.(['crest_',angleToken]) = crestFactor(:,is);
    metrics.(['kurtosis_',angleToken]) = kurtosisValue(:,is);
end
writetable(metrics,fullfile('results', ...
    ['paper_v3_three_channel_metrics',suffix,'.csv']));

showTime = time<=min(1,options.recordDuration);
colors = [0.08 0.08 0.08;0.05 0.35 0.85;0.85 0.15 0.10];
fig = figure('Visible','off','Color','w','Position',[50 40 1350 850]);
tiledlayout(nCase,nSensor,'TileSpacing','compact','Padding','compact');
for ic = 1:nCase
    for is = 1:nSensor
        nexttile;
        plot(time(showTime),responses{ic}.acceleration(showTime,is), ...
            'Color',colors(ic,:),'LineWidth',0.65);
        xlim([0,min(1,options.recordDuration)]); grid on; box off;
        xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
        title(sprintf('%s, sensor %.0f deg',caseNames{ic}, ...
            responses{ic}.sensorAnglesDeg(is)));
    end
end
sgtitle('Three-channel loaded 18+7 DOF effective transfer-path response');
exportgraphics(fig,fullfile('results', ...
    ['paper_v3_three_channel_time',suffix,'.png']),'Resolution',190);
close(fig);

fig = figure('Visible','off','Color','w','Position',[80 70 1180 850]);
tiledlayout(nSensor,1,'TileSpacing','compact','Padding','compact');
inMesh = frequency>=140 & frequency<=200;
for is = 1:nSensor
    nexttile;
    for ic = 1:nCase
        plot(frequency(inMesh),spectrum(inMesh,ic,is), ...
            'Color',colors(ic,:),'LineWidth',1.0); hold on;
    end
    for order = -3:3
        xline(p.kin.f_mesh+order*p.kin.f_c,':', ...
            'Color',[0.65 0.65 0.65],'HandleVisibility','off');
    end
    maximum = max(spectrum(inMesh,:,is),[],'all');
    xlim([140 200]); ylim([0,1.05*maximum]); grid on; box off;
    xlabel('Frequency (Hz)'); ylabel('Amplitude (m/s^2)');
    title(sprintf('Sensor %.0f deg: mesh band', ...
        responses{1}.sensorAnglesDeg(is)));
    legend(caseNames,'Location','best');
end
exportgraphics(fig,fullfile('results', ...
    ['paper_v3_three_channel_mesh_band',suffix,'.png']),'Resolution',190);
close(fig);

fprintf('Paper-v3 three-channel simulation complete.\n');
end
