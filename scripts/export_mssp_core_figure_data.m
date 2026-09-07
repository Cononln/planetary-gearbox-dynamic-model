function export_mssp_core_figure_data()
%EXPORT_MSSP_CORE_FIGURE_DATA Export existing MATLAB results for Python plots.
% This script performs data conversion only. It does not recompute experiments
% or create graphics.

rootDir = fileparts(fileparts(mfilename('fullpath')));
resultDir = fullfile(rootDir, 'results');
outDir = fullfile(resultDir, 'paper_figures_mssp', 'source_data');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Dynamic moving-path fields
s = load(fullfile(resultDir, 'loaded_path_am_pm.mat'), 'report');
f = s.report.pointField;
thetaDeg = rad2deg(f.phi_c(:));
amp = f.amplitude;
phaseDeg = rad2deg(f.phaseDifference);

Tpath = table(thetaDeg, ...
    squeeze(amp(1,:,1)).', squeeze(amp(1,:,2)).', phaseDeg(1,:).', ...
    squeeze(amp(2,:,1)).', squeeze(amp(2,:,2)).', phaseDeg(2,:).', ...
    f.reliability(1,:).', f.reliability(2,:).', ...
    'VariableNames', {'carrier_angle_deg', ...
    'amp_168_sensor0', 'amp_168_sensor90', 'phase_168_deg', ...
    'amp_1848_sensor0', 'amp_1848_sensor90', 'phase_1848_deg', ...
    'reliability_168', 'reliability_1848'});
writetable(Tpath, fullfile(outDir, 'dynamic_path_field.csv'));

%% Encoder-hidden validation trajectories and path truth
s = load(fullfile(resultDir, 'real_tacholess_joint_phase.mat'), 'report');
r = s.report;
stride = max(1, round(r.fs / 100));
idx = (1:stride:numel(r.t)).';
Tfreq = table(r.t(idx), r.truth.meshFrequencyHz(idx), ...
    r.joint.instantaneousMeshFrequencyHz(idx), ...
    'VariableNames', {'time_s', 'encoder_mesh_frequency_hz', ...
    'blind_mesh_frequency_hz'});
writetable(Tfreq, fullfile(outDir, 'encoder_hidden_frequency.csv'));

nb = 180;
edges = linspace(0, 2*pi, nb + 1);
centres = (edges(1:end-1) + edges(2:end)) / 2;
theta = mod(r.truth.carrierPhase, 2*pi);
oracleCross = r.oraclePathFull(:,2,1) .* conj(r.oraclePathFull(:,1,1));
blindCross = r.hybrid.pathPhasor(:,2,1) .* conj(r.hybrid.pathPhasor(:,1,1));
oraclePhase = nan(nb,1);
blindPhase = nan(nb,1);
support = zeros(nb,1);
for b = 1:nb
    mask = theta >= edges(b) & theta < edges(b+1);
    support(b) = nnz(mask);
    if support(b) > 0
        oraclePhase(b) = angle(sum(oracleCross(mask) ./ max(abs(oracleCross(mask)), eps)));
        blindPhase(b) = angle(sum(blindCross(mask) ./ max(abs(blindCross(mask)), eps)));
    end
end
Tphase = table(rad2deg(centres(:)), rad2deg(oraclePhase), ...
    rad2deg(blindPhase), support, ...
    'VariableNames', {'carrier_angle_deg', 'encoder_oracle_path_deg', ...
    'blind_path_deg', 'sample_count'});
writetable(Tphase, fullfile(outDir, 'encoder_hidden_path.csv'));
writetable(r.metrics, fullfile(outDir, 'encoder_hidden_metrics.csv'));

%% Sun-fault diagnostic source data
s = load(fullfile(resultDir, 'real_phase_aligned_diagnosis.mat'), 'diagnosis');
d = s.diagnosis;
keep = d.frequency <= 105;
Tspec = table(d.frequency(keep), d.spectrum(keep,1), d.spectrum(keep,2), ...
    d.spectrum(keep,3), 'VariableNames', {'frequency_hz', 'single_sensor', ...
    'naive_fusion', 'phase_aligned_fusion'});
writetable(Tspec, fullfile(outDir, 'sun_fault_envelope_spectrum.csv'));
writetable(d.metrics, fullfile(outDir, 'sun_fault_metrics.csv'));

angleDeg = linspace(0, 360, size(d.cycleMatrix{1},2) + 1);
angleDeg(end) = [];
for k = 1:3
    cycleTable = array2table(d.cycleMatrix{k});
    cycleTable = addvars(cycleTable, (1:size(d.cycleMatrix{k},1)).', ...
        'Before', 1, 'NewVariableNames', 'cycle_index');
    writetable(cycleTable, fullfile(outDir, sprintf('sun_fault_cycle_matrix_%d.csv', k)));
end
writematrix(angleDeg(:), fullfile(outDir, 'sun_fault_cycle_angle_deg.csv'));

pattern = [d.tpsvdPattern{1}(:), d.tpsvdPattern{2}(:), d.tpsvdPattern{3}(:)];
Tpattern = array2table(pattern, 'VariableNames', ...
    {'single_sensor', 'naive_fusion', 'phase_aligned_fusion'});
writetable(Tpattern, fullfile(outDir, 'sun_fault_tpsvd_pattern.csv'));

maxRank = max(cellfun(@numel, d.singularValue));
sv = nan(maxRank,3);
for k = 1:3
    sv(1:numel(d.singularValue{k}),k) = d.singularValue{k}(:);
end
Tsv = array2table([(1:maxRank).', sv], 'VariableNames', ...
    {'rank', 'single_sensor', 'naive_fusion', 'phase_aligned_fusion'});
writetable(Tsv, fullfile(outDir, 'sun_fault_singular_values.csv'));

fprintf('Exported MSSP core figure source data to %s\n', outDir);
end
