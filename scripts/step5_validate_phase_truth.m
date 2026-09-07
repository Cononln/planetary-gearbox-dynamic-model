function report = step5_validate_phase_truth
%STEP5_VALIDATE_PHASE_TRUTH Generate and validate phase/group-delay truth.

setup_paths;
p = pg_parameters;
model = pg_ring_modal_model(p);
frequencyHz = (1600:10:5200).';
phi = linspace(0,2*pi,181);
truth = pg_phase_truth_field(frequencyHz,phi,p,model,1);

assert(isequal(size(truth.H0),[numel(frequencyHz),numel(phi)]));
assert(all(isfinite(truth.H0(:))) && all(isfinite(truth.H90(:))));
assert(nnz(truth.validMask)>0.5*numel(truth.validMask));

periodicH0 = norm(truth.H0(:,end)-truth.H0(:,1))/norm(truth.H0(:,1));
periodicH90 = norm(truth.H90(:,end)-truth.H90(:,1))/norm(truth.H90(:,1));
assert(max(periodicH0,periodicH90)<1e-10);

% Construct a known broadband excitation map. Since the same excitation drives
% both paths, their cross-phase is exactly the model truth.
excitation = exp(-0.5*((frequencyHz-1840)/90).^2) + ...
             0.8*exp(-0.5*((frequencyHz-3030)/110).^2) + ...
             0.6*exp(-0.5*((frequencyHz-4980)/80).^2);
excitation = repmat(excitation,1,numel(phi));
S0 = truth.H0.*excitation;
S90 = truth.H90.*excitation;
S90Aligned = pg_apply_phase_truth(S90,truth);

before = angle(S90.*conj(S0));
after = angle(S90Aligned.*conj(S0));
weight = truth.reliability.*abs(S0).*abs(S90);
use = truth.validMask & weight>0;
rmsBefore = circularRms(before(use),weight(use));
rmsAfter = circularRms(after(use),weight(use));
assert(rmsAfter<1e-10);
assert(rmsBefore>deg2rad(5));

% Strong static comparator: use the complete frequency correction at phi=0
% for every angle. It still cannot remove a carrier-angle-dependent path.
staticCorrection = repmat(truth.phaseDifferenceWrapped(:,1),1,numel(phi));
S90Static = S90.*exp(-1i*staticCorrection);
staticResidual = angle(S90Static.*conj(S0));
rmsStatic = circularRms(staticResidual(use),weight(use));
assert(rmsStatic>100*rmsAfter+deg2rad(1));

bandSpanUs = zeros(1,size(truth.bandsHz,1));
for ib = 1:size(truth.bandsHz,1)
    values = 1e6*truth.bandDelaySec(ib,:);
    bandSpanUs(ib) = max(values,[],'omitnan')-min(values,[],'omitnan');
end
assert(any(bandSpanUs>1));

report.frequencyRangeHz = [frequencyHz(1),frequencyHz(end)];
report.carrierAngleStepDeg = rad2deg(phi(2)-phi(1));
report.validFraction = nnz(truth.validMask)/numel(truth.validMask);
report.periodicResidual = max(periodicH0,periodicH90);
report.rmsPhaseBeforeDeg = rad2deg(rmsBefore);
report.rmsPhaseStaticDeg = rad2deg(rmsStatic);
report.rmsPhaseTruthDeg = rad2deg(rmsAfter);
report.bandDelaySpanUs = bandSpanUs;
report.bandsHz = truth.bandsHz;

save(fullfile('results','phase_truth_field.mat'),'truth','-v7.3');
delayTable = table(rad2deg(phi(:)), ...
    1e6*truth.bandDelaySec(1,:).',1e6*truth.bandDelaySec(2,:).', ...
    1e6*truth.bandDelaySec(3,:).', ...
    'VariableNames',{'carrier_angle_deg','delay_1830band_us', ...
    'delay_3030band_us','delay_4980band_us'});
writetable(delayTable,fullfile('results','phase_truth_band_delay.csv'));

groupDelayUs = 1e6*truth.groupDelaySec;
groupDelayUs(~truth.validMask) = NaN;
finiteDelay = groupDelayUs(isfinite(groupDelayUs));
sortedAbsDelay = sort(abs(finiteDelay));
clipLim = sortedAbsDelay(max(1,ceil(0.98*numel(sortedAbsDelay))));
groupDelayUs = max(min(groupDelayUs,clipLim),-clipLim);

fig = figure('Visible','off','Color','w','Position',[100 100 1040 850]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile;
imagesc(rad2deg(phi),frequencyHz,rad2deg(truth.phaseDifferenceWrapped));
axis xy; colorbar; caxis([-180 180]);
xlabel('Carrier angle (deg)'); ylabel('Frequency (Hz)');
title('True dual-sensor phase difference (deg)');
nexttile;
imagesc(rad2deg(phi),frequencyHz,groupDelayUs);
axis xy; colorbar; caxis([-clipLim clipLim]);
xlabel('Carrier angle (deg)'); ylabel('Frequency (Hz)');
title(sprintf('True group delay (microseconds, clipped at 98%% = %.1f)',clipLim));
nexttile;
plot(rad2deg(phi),1e6*truth.bandDelaySec,'LineWidth',1.1);
xlim([0 360]); xlabel('Carrier angle (deg)');
ylabel('Band-weighted group delay (microseconds)');
legend('1.75-1.92 kHz','2.92-3.14 kHz','4.89-5.06 kHz', ...
    'Location','best'); grid on;
title('Periodic angle-dependent delay fields');
exportgraphics(fig,fullfile('results','step5_phase_truth.png'),'Resolution',180);
close(fig);

fprintf('Stage 5 validation PASSED\n');
fprintf('  valid field fraction       : %.3f\n',report.validFraction);
fprintf('  2pi periodic residual      : %.3e\n',report.periodicResidual);
fprintf('  phase RMS before/static/truth: %.2f / %.2f / %.3e deg\n', ...
    report.rmsPhaseBeforeDeg,report.rmsPhaseStaticDeg,report.rmsPhaseTruthDeg);
fprintf('  band delay spans (us)      :'); fprintf(' %.2f',bandSpanUs); fprintf('\n');
end

function value = circularRms(angleValue,weight)
weight = weight(:);
angleValue = mod(angleValue(:)+pi,2*pi)-pi;
value = sqrt(sum(weight.*angleValue.^2)/sum(weight));
end
