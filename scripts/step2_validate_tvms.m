function report = step2_validate_tvms
%STEP2_VALIDATE_TVMS Validate healthy, 25%, and 50% broken-tooth TVMS.

setup_paths;
p = pg_parameters;
healthy = pg_fault_case('none',0,p);
sun25 = pg_fault_case('sun',0.25,p);
sun50 = pg_fault_case('sun',0.50,p);
planet25 = pg_fault_case('planet',0.25,p);
planet50 = pg_fault_case('planet',0.50,p);

% One mesh period is sufficient to verify healthy mean calibration.
nCycle = 8192;
tCycle = (0:nCycle-1)'/(nCycle*p.kin.f_mesh);
kHealthyCycle = pg_tvms(tCycle,p,healthy);
meanSp = mean(kHealthyCycle.sunPlanet,1);
meanRp = mean(kHealthyCycle.ringPlanet,1);
meanErrorSp = max(abs(meanSp-p.mesh.sunPlanet.kMean))/p.mesh.sunPlanet.kMean;
meanErrorRp = max(abs(meanRp-p.mesh.ringPlanet.kMean))/p.mesh.ringPlanet.kMean;
assert(meanErrorSp < 5e-4);
assert(meanErrorRp < 5e-4);
assert(all(kHealthyCycle.sunPlanet(:)>0));
assert(all(kHealthyCycle.ringPlanet(:)>0));

% A quarter second contains six sun-fault events and multiple planet events.
t = (0:1/p.model.fs:0.25-1/p.model.fs)';
kH = pg_tvms(t,p,healthy);
kS25 = pg_tvms(t,p,sun25);
kS50 = pg_tvms(t,p,sun50);
kP25 = pg_tvms(t,p,planet25);
kP50 = pg_tvms(t,p,planet50);

sumHsp = sum(kH.sunPlanet,2);
sumS25 = sum(kS25.sunPlanet,2);
sumS50 = sum(kS50.sunPlanet,2);
dropS25 = sumHsp-sumS25;
dropS50 = sumHsp-sumS50;
assert(max(dropS25)>0 && max(dropS50)>max(dropS25));
assert(abs(max(dropS50)/max(dropS25)-2)<0.03);

% The total mesh is unchanged outside the damaged tooth's contact window,
% so a global constant multiplier is explicitly rejected.
ratioS25 = sumS25./sumHsp;
assert(max(ratioS25)>0.999999);
assert(min(ratioS25)<0.999);
assert(std(ratioS25)>1e-3);

% Planet 1 tooth 1 must affect both its sun and ring mesh; planets 2/3 stay healthy.
assert(any(kP25.sunPlanet(:,1)<kH.sunPlanet(:,1)));
assert(any(kP25.ringPlanet(:,1)<kH.ringPlanet(:,1)));
assert(max(abs(kP25.sunPlanet(:,2:3)-kH.sunPlanet(:,2:3)),[],'all')<1e-9);
assert(max(abs(kP25.ringPlanet(:,2:3)-kH.ringPlanet(:,2:3)),[],'all')<1e-9);
planetDrop25 = sum(kH.sunPlanet+kH.ringPlanet,2)- ...
               sum(kP25.sunPlanet+kP25.ringPlanet,2);
planetDrop50 = sum(kH.sunPlanet+kH.ringPlanet,2)- ...
               sum(kP50.sunPlanet+kP50.ringPlanet,2);
assert(max(planetDrop50)>max(planetDrop25));

% The aggregate sun-fault pattern repeats after seven mesh cycles = 1/24 s.
tProbe = linspace(0,1/p.kin.f_sun_fault,1500).';
kA = pg_tvms(tProbe,p,sun25);
kB = pg_tvms(tProbe+1/p.kin.f_sun_fault,p,sun25);
sunPeriodResidual = norm(sum(kA.sunPlanet,2)-sum(kB.sunPlanet,2))/ ...
                    norm(sum(kA.sunPlanet,2));
assert(sunPeriodResidual<1e-12);

% Full matrices remain positive definite under the most severe frozen case.
sampleTimes = linspace(0,1/p.kin.f_sun_fault,12);
minK = inf;
for it = 1:numel(sampleTimes)
    [~,~,K] = pg_assemble_system(sampleTimes(it),p,sun50);
    minK = min(minK,min(eig(K)));
end
assert(minK>0);

report.meanErrorSunPlanet = meanErrorSp;
report.meanErrorRingPlanet = meanErrorRp;
report.sun25MissingWidthMm = 1e3*sun25.missingFaceWidth;
report.sun50MissingWidthMm = 1e3*sun50.missingFaceWidth;
report.radialDepthMm = 1e3*sun25.radialRemovalDepth;
report.maxSunDrop25 = max(dropS25);
report.maxSunDrop50 = max(dropS50);
report.sunPeriodResidual = sunPeriodResidual;
report.minSystemStiffnessEigenvalue = minK;

fig = figure('Visible','off','Color','w','Position',[100 100 920 620]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;
plot(t, sumHsp/3/p.mesh.sunPlanet.kMean,'k','LineWidth',1.1); hold on;
plot(t, sumS25/3/p.mesh.sunPlanet.kMean,'b','LineWidth',1.0);
plot(t, sumS50/3/p.mesh.sunPlanet.kMean,'r','LineWidth',1.0);
xlabel('Time (s)'); ylabel('Normalized SP stiffness');
legend('Healthy','25% broken width','50% broken width','Location','best');
title('Aggregate sun-planet TVMS'); grid on;
nexttile;
totalHP1 = kH.sunPlanet(:,1)+kH.ringPlanet(:,1);
totalP25 = kP25.sunPlanet(:,1)+kP25.ringPlanet(:,1);
totalP50 = kP50.sunPlanet(:,1)+kP50.ringPlanet(:,1);
plot(t,totalHP1/mean(totalHP1),'k','LineWidth',1.1); hold on;
plot(t,totalP25/mean(totalHP1),'b','LineWidth',1.0);
plot(t,totalP50/mean(totalHP1),'r','LineWidth',1.0);
xlabel('Time (s)'); ylabel('Normalized P1 mesh stiffness');
legend('Healthy','25% broken width','50% broken width','Location','best');
title('Planet 1 tooth fault in its two meshes'); grid on;
exportgraphics(fig,fullfile('results','step2_tvms.png'),'Resolution',180);
close(fig);

fprintf('Stage 2 validation PASSED\n');
fprintf('  healthy mean errors SP/RP : %.3e / %.3e\n',meanErrorSp,meanErrorRp);
fprintf('  missing widths 25/50      : %.2f / %.2f mm\n', ...
    report.sun25MissingWidthMm,report.sun50MissingWidthMm);
fprintf('  radial removal depth      : %.3f mm (to root circle)\n', ...
    report.radialDepthMm);
fprintf('  aggregate sun period      : %.6f s (%.3f Hz)\n', ...
    1/p.kin.f_sun_fault,p.kin.f_sun_fault);
end
