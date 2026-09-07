function report = validate_root_crack_model
%VALIDATE_ROOT_CRACK_MODEL Validate Shen-style sun-root crack TVMS.

setup_paths;
p = pg_parameters;
healthy = pg_fault_case('none',0,p);
paperDepthsMm = [0.5,1.5,2.5];
paperModuleMm = 2.0;
depthsMm = paperDepthsMm*(1e3*p.gear.module/paperModuleMm);
faceCoverage = 0.5;
faults = arrayfun(@(q) pg_root_crack_case(q,p, ...
    'FaceCoverage',faceCoverage),depthsMm, ...
    'UniformOutput',false);

for k = 1:numel(faults)-1
    assert(min(faults{k+1}.energyLookup.stiffnessRatio) < ...
        min(faults{k}.energyLookup.stiffnessRatio));
end

duration = 1.0;
t = (0:1/p.model.fs:duration-1/p.model.fs).';
kHealthy = pg_tvms(t,p,healthy);
aggregateHealthy = sum(kHealthy.sunPlanet,2);
aggregate = zeros(numel(t),numel(faults));
deficit = zeros(size(aggregate));
for k = 1:numel(faults)
    tvms = pg_tvms(t,p,faults{k});
    aggregate(:,k) = sum(tvms.sunPlanet,2);
    deficit(:,k) = aggregateHealthy-aggregate(:,k);
    assert(max(deficit(:,k))>0);
end
assert(all(max(deficit(:,2:end),[],1)>max(deficit(:,1:end-1),[],1)));

% Fault pattern must repeat at the three-planet sun-fault frequency.
tProbe = linspace(0,1/p.kin.f_sun_fault,1600).';
kA = pg_tvms(tProbe,p,faults{2});
kB = pg_tvms(tProbe+1/p.kin.f_sun_fault,p,faults{2});
periodResidual = norm(sum(kA.sunPlanet,2)-sum(kB.sunPlanet,2))/ ...
    norm(sum(kA.sunPlanet,2));
assert(periodResidual<1e-11);

n = numel(t);
window = 0.5-0.5*cos(2*pi*(0:n-1)'/(n-1));
frequency = (0:floor(n/2))'*p.model.fs/n;
deficitSpectrum = zeros(numel(frequency),numel(faults));
for k = 1:numel(faults)
    X = fft((deficit(:,k)-mean(deficit(:,k))).*window);
    deficitSpectrum(:,k) = 2*abs(X(1:numel(frequency)))/sum(window);
end

targetFrequencies = [p.kin.f_sun_fault,p.kin.f_mesh-p.kin.f_sun_fault, ...
    p.kin.f_mesh,p.kin.f_mesh+p.kin.f_sun_fault];
targetAmplitude = zeros(numel(faults),numel(targetFrequencies));
for j = 1:numel(targetFrequencies)
    [~,idx] = min(abs(frequency-targetFrequencies(j)));
    targetAmplitude(:,j) = deficitSpectrum(idx,:).';
end

report.depthMm = depthsMm;
report.paperEquivalentDepthMm = paperDepthsMm;
report.faceCoverage = faceCoverage;
report.depthToModuleRatio = cellfun(@(f) f.depthToModuleRatio,faults);
report.minimumPairRatio = cellfun(@(f) ...
    min(f.energyLookup.stiffnessRatio),faults);
report.maximumAggregateDrop = max(deficit,[],1);
report.periodResidual = periodResidual;
report.faultFrequencyHz = p.kin.f_sun_fault;
report.targetFrequenciesHz = targetFrequencies;
report.targetAmplitudeNm = targetAmplitude;

if ~exist('results','dir'); mkdir('results'); end
metrics = table(paperDepthsMm.',depthsMm.', ...
    repmat(faceCoverage,numel(depthsMm),1),report.depthToModuleRatio.', ...
    report.minimumPairRatio.', ...
    report.maximumAggregateDrop.',targetAmplitude(:,1), ...
    targetAmplitude(:,2),targetAmplitude(:,3),targetAmplitude(:,4), ...
    'VariableNames',{'paper_crack_depth_mm','scaled_crack_depth_mm', ...
    'face_coverage','depth_to_module_ratio', ...
    'minimum_pair_stiffness_ratio', ...
    'maximum_aggregate_drop_Npm','fault_order_amplitude_Npm', ...
    'lower_sideband_amplitude_Npm','mesh_amplitude_Npm', ...
    'upper_sideband_amplitude_Npm'});
writetable(metrics,fullfile('results','root_crack_tvms_metrics.csv'));

colors = [0.20 0.45 0.75;0.90 0.55 0.10;0.75 0.15 0.18];
fig = figure('Visible','off','Color','w','Position',[80 60 1120 760]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
for k = 1:numel(faults)
    plot(100*faults{k}.energyLookup.contactFraction, ...
        faults{k}.energyLookup.stiffnessRatio,'Color',colors(k,:), ...
        'LineWidth',1.4); hold on;
end
xlabel('Contact progress (%)'); ylabel('Cracked/healthy pair stiffness');
title('(a) Energy-based local stiffness ratio');
legend(compose('q = %.3f mm',depthsMm),'Location','best'); grid on;

nexttile;
plot(t,aggregateHealthy/(3*p.mesh.sunPlanet.kMean),'k','LineWidth',0.9); hold on;
for k = 1:numel(faults)
    plot(t,aggregate(:,k)/(3*p.mesh.sunPlanet.kMean), ...
        'Color',colors(k,:),'LineWidth',0.8);
end
xlim([0 0.12]); xlabel('Time (s)');
ylabel('Normalized aggregate SP stiffness');
title('(b) Tooth-indexed sun-planet TVMS'); grid on;

nexttile;
plot(t,deficit(:,1)/1e6,'Color',colors(1,:),'LineWidth',0.8); hold on;
plot(t,deficit(:,2)/1e6,'Color',colors(2,:),'LineWidth',0.8);
plot(t,deficit(:,3)/1e6,'Color',colors(3,:),'LineWidth',0.8);
xlim([0 0.12]); xlabel('Time (s)'); ylabel('Stiffness deficit (MN/m)');
title(sprintf('(c) Repetition at %.1f Hz',p.kin.f_sun_fault)); grid on;

nexttile;
inBand = frequency<=240;
for k = 1:numel(faults)
    plot(frequency(inBand),deficitSpectrum(inBand,k)/1e6, ...
        'Color',colors(k,:),'LineWidth',1.0); hold on;
end
xline(p.kin.f_sun_fault,':','Color',[0.35 0.35 0.35]);
xline(p.kin.f_mesh,'--','Color',[0.35 0.35 0.35]);
xlabel('Frequency (Hz)'); ylabel('Deficit amplitude (MN/m)');
title('(d) Fault-order and mesh-sideband content'); grid on;
exportgraphics(fig,fullfile('results','root_crack_tvms_validation.png'), ...
    'Resolution',200);
close(fig);

save(fullfile('results','root_crack_tvms_validation.mat'), ...
    'report','faults','-v7.3');
fprintf('Root-crack TVMS validation PASSED\n');
fprintf('  scaled q (mm)           : %.3f / %.3f / %.3f\n',depthsMm);
fprintf('  crack face coverage     : %.1f percent\n',100*faceCoverage);
fprintf('  minimum pair ratios     : %.4f / %.4f / %.4f\n', ...
    report.minimumPairRatio);
fprintf('  repetition frequency    : %.3f Hz\n',p.kin.f_sun_fault);
fprintf('  period residual         : %.3e\n',periodResidual);
end
