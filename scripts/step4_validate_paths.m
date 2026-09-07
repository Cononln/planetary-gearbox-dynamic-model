function report = step4_validate_paths
%STEP4_VALIDATE_PATHS Validate 24-DOF coupling and moving-source paths.

setup_paths;
p = pg_parameters;
model = pg_ring_modal_model(p);
healthy = pg_fault_case('none',0,p);

angles = [0,pi/9,pi/2,2*pi];
couplingNorm = zeros(size(angles));
for k = 1:numel(angles)
    [M,C,K,meta] = pg_coupled_matrices(angles(k),p,model);
    assert(isequal(size(M),[24,24]));
    assert(norm(M-M.','fro')<1e-12);
    assert(norm(C-C.','fro')<1e-10);
    assert(norm(K-K.','fro')<1e-6);
    assert(min(eig(K))>0);
    couplingNorm(k) = norm(K(1:18,19:24),'fro');
    assert(couplingNorm(k)>0);
    for i = 1:p.model.nPlanet
        assert(numel(meta.coupling(i).g)==24);
    end
end

[M0,C0,K0] = pg_coupled_matrices(0,p,model);
[M2,C2,K2] = pg_coupled_matrices(2*pi,p,model);
periodicResidual = max([norm(M2-M0,'fro')/norm(M0,'fro'), ...
    norm(C2-C0,'fro')/norm(C0,'fro'), ...
    norm(K2-K0,'fro')/norm(K0,'fro')]);
assert(periodicResidual<1e-12);

Cacc = pg_sensor_observation(p,model);
assert(isequal(size(Cacc),[2,24]));
assert(norm(Cacc(1,p.map.ring(1:2))-[1,0])<1e-12);
assert(norm(Cacc(2,p.map.ring(1:2))-[0,1])<1e-12);

% Check the time-domain interface without integrating a long experiment.
dx0 = pg_rhs_coupled(0,zeros(48,1),p,healthy,model);
assert(isequal(size(dx0),[48,1]) && all(isfinite(dx0)));

phi = linspace(0,2*pi,721);
frequencyHz = 1840;
H0 = zeros(size(phi));
H90 = zeros(size(phi));
for k = 1:numel(phi)
    path = pg_path_frf(frequencyHz,phi(k),p,model);
    H0(k) = path.H(1,1);
    H90(k) = path.H(2,1);
end
phaseDifference = unwrap(angle(H90))-unwrap(angle(H0));
valid = abs(H0)>max(abs(H0))*1e-5 & abs(H90)>max(abs(H90))*1e-5;
phaseSpan = range(phaseDifference(valid));
amplitudeCv0 = std(abs(H0))/mean(abs(H0));
amplitudeCv90 = std(abs(H90))/mean(abs(H90));
assert(phaseSpan>deg2rad(10));
assert(amplitudeCv0>0.05 && amplitudeCv90>0.05);

report.coupledDof = 24;
report.periodicResidual = periodicResidual;
report.couplingBlockNormRange = [min(couplingNorm),max(couplingNorm)];
report.pathFrequencyHz = frequencyHz;
report.phaseSpanDeg = rad2deg(phaseSpan);
report.amplitudeCv = [amplitudeCv0,amplitudeCv90];

fig = figure('Visible','off','Color','w','Position',[100 100 960 680]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;
plot(rad2deg(phi),20*log10(abs(H0)/max(abs(H0))),'b','LineWidth',1.1); hold on;
plot(rad2deg(phi),20*log10(abs(H90)/max(abs(H90))),'r','LineWidth',1.1);
xlim([0 360]); xlabel('Carrier angle (deg)');
ylabel('Relative path magnitude (dB)');
legend('0 deg sensor','90 deg sensor','Location','best'); grid on;
title(sprintf('Moving mesh 1 path at %.0f Hz',frequencyHz));
nexttile;
plot(rad2deg(phi(valid)),rad2deg(wrapToPiLocal(phaseDifference(valid))), ...
    'k','LineWidth',1.1);
xlim([0 360]); ylim([-180 180]);
xlabel('Carrier angle (deg)'); ylabel('H90-H0 phase (deg)');
title('Dual-sensor physical path phase difference'); grid on;
exportgraphics(fig,fullfile('results','step4_paths.png'),'Resolution',180);
close(fig);

fprintf('Stage 4 validation PASSED\n');
fprintf('  coupled DOF             : %d\n',report.coupledDof);
fprintf('  2pi periodic residual   : %.3e\n',periodicResidual);
fprintf('  phase span at %.0f Hz    : %.2f deg\n',frequencyHz,report.phaseSpanDeg);
fprintf('  path magnitude CV 0/90  : %.3f / %.3f\n',amplitudeCv0,amplitudeCv90);
end

function y = wrapToPiLocal(x)
y = mod(x+pi,2*pi)-pi;
end
