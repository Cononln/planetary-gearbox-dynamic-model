function export_ch3_path_truth_csv(inputMat, outputCsv)
%EXPORT_CH3_PATH_TRUTH_CSV Export legacy 18+7 path-truth amplitudes/phases.
% This is a read-only post-processing helper. It does not run the FE path.
arguments
    inputMat (1,1) string
    outputCsv (1,1) string
end
S = load(inputMat, 'pathTruth', 'eventBook');
T = S.pathTruth;
f = T.frequencyHz(:);
target = [168, 336, 504, 1848];
[~, idx] = arrayfun(@(x) min(abs(f-x)), target);
idx = unique(idx, 'stable');
eventAngle = mod(rad2deg(T.carrierAngleRad(:)), 360);
eventId = T.eventId(:);
planet = T.sourcePlanetIndex(:);
nEvent = numel(eventId);
nSensor = size(T.amplitude, 3);
labels = string(S.eventBook.sensorLabels(:));
N = numel(idx)*nEvent*nSensor;
frequencyOut = zeros(N,1);
eventOut = zeros(N,1);
angleOut = zeros(N,1);
planetOut = zeros(N,1);
sensorOut = strings(N,1);
amplitudeOut = zeros(N,1);
relativeOut = zeros(N,1);
phaseOut = zeros(N,1);
k = 0;
for jf = 1:numel(idx)
    fi = idx(jf);
    for ie = 1:nEvent
        for is = 1:nSensor
            k = k + 1;
            frequencyOut(k) = f(fi);
            eventOut(k) = eventId(ie);
            angleOut(k) = eventAngle(ie);
            planetOut(k) = planet(ie);
            sensorOut(k) = labels(is);
            amplitudeOut(k) = T.amplitude(fi,ie,is);
            relativeOut(k) = T.relativeAmplitude(fi,ie,is);
            phaseOut(k) = T.phaseDifferenceWrapped(fi,ie,is);
        end
    end
end
out = table(frequencyOut,eventOut,angleOut,planetOut,sensorOut, ...
    amplitudeOut,relativeOut,phaseOut, 'VariableNames', ...
    {'frequency_Hz','event_id','carrier_angle_deg','source_planet', ...
     'sensor_label','amplitude','relative_amplitude','phase_difference_deg'});
writetable(out, outputCsv);
fprintf('Exported %d rows to %s\n', height(out), outputCsv);
end
