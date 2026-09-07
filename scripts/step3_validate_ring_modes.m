function report = step3_validate_ring_modes
%STEP3_VALIDATE_RING_MODES Validate frequency calibration and analytical basis.

setup_paths;
p = pg_parameters;
model = pg_ring_modal_model(p);
expectedHz = [1832.565469,1845.781068,3010.794375, ...
              3053.891626,4976.103825,4976.220685];
assert(max(abs(model.frequencyHz-expectedHz))<1e-9);
assert(isequal(model.waveNumber,[2,2,3,3,4,4]));

theta = (0:4095)'*(2*pi/4096);
Psi = pg_ring_mode_shape(theta,model);
gram = (Psi.'*Psi)/numel(theta);
orthogonalityResidual = norm(gram-0.5*eye(6),'fro');
assert(orthogonalityResidual<1e-12);

assert(norm(model.M-model.M.','fro')==0);
assert(norm(model.C-model.C.','fro')==0);
assert(norm(model.K-model.K.','fro')==0);
assert(all(diag(model.M)>0) && all(diag(model.K)>0));

% Confirm the expected three moving input columns and point observations.
[B0,thetaMesh] = pg_ring_input_matrix(0,p,model);
C0 = pg_ring_observation(0,model);
C90 = pg_ring_observation(pi/2,model);
assert(isequal(size(B0),[6,3]));
assert(isequal(size(C0),[1,6]) && isequal(size(C90),[1,6]));
assert(max(abs(diff(sort(mod(thetaMesh,2*pi)))-2*pi/3))<1e-12);

% At the 0-degree collocated point the sine members have nodes; the response
% is nevertheless finite and displays the observable cosine-member resonances.
f = linspace(500,5500,12000);
hPoint = zeros(size(f));
bPoint = model.radialMeshProjection*pg_ring_mode_shape(0,model).';
for k = 1:numel(f)
    D = pg_ring_dynamic_stiffness(2*pi*f(k),model);
    hPoint(k) = C0*(D\bPoint);
end
frfFinite = all(isfinite(hPoint));
assert(frfFinite);

pairSplitHz = [diff(model.frequencyHz(1:2)), ...
               diff(model.frequencyHz(3:4)), ...
               diff(model.frequencyHz(5:6))];
assert(all(pairSplitHz>0));

report.sourceModes = model.sourceModes;
report.frequencyHz = model.frequencyHz;
report.waveNumber = model.waveNumber;
report.pairSplitHz = pairSplitHz;
report.orthogonalityResidual = orthogonalityResidual;
report.modalMassKg = model.modalMass;

fig = figure('Visible','off','Color','w','Position',[100 100 980 650]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;
plot(rad2deg(theta),Psi,'LineWidth',1.0);
xlim([0 360]); xticks(0:60:360);
xlabel('Circumferential angle (deg)'); ylabel('Radial modal amplitude');
title('Frequency-calibrated analytical ring basis');
legend('7: cos2','8: sin2','9: cos3','10: sin3','11: cos4','12: sin4', ...
       'NumColumns',3,'Location','best'); grid on;
nexttile;
semilogy(f,abs(hPoint),'k','LineWidth',1.1); hold on;
for r = 1:numel(model.frequencyHz)
    xline(model.frequencyHz(r),'--','Color',[0.75 0.15 0.15]);
end
xlabel('Frequency (Hz)'); ylabel('|H_{rr}| (m/N)');
title('Radial point receptance at 0 degrees'); grid on;
exportgraphics(fig,fullfile('results','step3_ring_modes.png'),'Resolution',180);
close(fig);

fprintf('Stage 3 validation PASSED\n');
fprintf('  source modes       :'); fprintf(' %d',model.sourceModes); fprintf('\n');
fprintf('  frequencies (Hz)   :'); fprintf(' %.3f',model.frequencyHz); fprintf('\n');
fprintf('  pair splits (Hz)   :'); fprintf(' %.3f',pairSplitHz); fprintf('\n');
fprintf('  basis Gram residual: %.3e\n',orthogonalityResidual);
end
