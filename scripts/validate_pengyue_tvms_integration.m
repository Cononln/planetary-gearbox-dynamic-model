function report = validate_pengyue_tvms_integration(options)
%VALIDATE_PENGYUE_TVMS_INTEGRATION Validate TVMS before dynamics.
% Generates healthy and q=3 mm, gamma=45 deg SP TVMS only.  The output
% directory is configurable because some installations keep the repository
% read-only.

arguments
    options.OutputDir (1,1) string = ""
    options.Samples (1,1) double {mustBeInteger,mustBeGreaterThan(options.Samples,1000)} = 6001
end

setup_paths;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
if strlength(options.OutputDir)==0
    outDir = fullfile(projectRoot,'results','pengyue_tvms_validation');
else
    outDir = char(options.OutputDir);
end
if ~isfolder(outDir); mkdir(outDir); end
p = pg_parameters_pengyue_fig314;
fHealthy = struct('type','none','tooth',1,'crackDepth',0,'crackAngle',0);
fCrack = struct('type','sun','tooth',1,'crackDepth',3e-3, ...
    'crackAngle',deg2rad(45),'crackAngleDeg',45);

[pH,modelH] = pg_prepare_pengyue_tvms(p,fHealthy);
[pF,modelF] = pg_prepare_pengyue_tvms(p,fCrack);
u = linspace(0,3,options.Samples).';
t = u/p.kin.f_mesh;
[kH,dH] = pg_eval_pengyue_sp_tvms(t,pH,fHealthy,modelH);
[kF,dF] = pg_eval_pengyue_sp_tvms(t,pF,fCrack,modelF);

% Verify the actual assembly sees the same 18-DOF structure and TVMS values.
[M,C,K,meta] = pg_assemble_pengyue_system(t(1),pF,fCrack); %#ok<ASGLU>
assert(isequal(size(M),[18 18]),'Expected an 18x18 mass matrix.');
assert(max(abs(meta.meshStiffness.sunPlanet(:)-kF(1,:).')) < 1e-6, ...
    'Prepared TVMS did not reach the assembly.');

activeH = dH.activePairCount(:,1);
doubleFraction = mean(activeH==2);
lossFraction = 100*(1-min(kF(:,1)./max(kH(:,1),eps)));
fprintf('Peng Yue SP TVMS validation PASSED\n');
fprintf('  DOF = %d\n',size(M,1));
fprintf('  fm = %.12f Hz, epsilon_alpha = %.9f\n',p.kin.f_mesh,modelH.contactRatio);
fprintf('  healthy k_sp1: %.9g .. %.9g N/m (mean %.9g)\n', ...
    min(kH(:,1)),max(kH(:,1)),mean(kH(:,1)));
fprintf('  crack q=3mm/45deg k_sp1: %.9g .. %.9g N/m (mean %.9g)\n', ...
    min(kF(:,1)),max(kF(:,1)),mean(kF(:,1)));
fprintf('  double-pair fraction = %.6f; crack loss metric = %.3f %%\n', ...
    doubleFraction,lossFraction);

T = table(u,t,kH(:,1),kF(:,1),activeH,dF.faultPairCount(:,1), ...
    dF.faultStiffnessLoss(:,1), ...
    'VariableNames',{'mesh_cycle','time_s','healthy_Npm','crack_q3_gamma45_Npm', ...
    'active_pair_count','fault_pair_count','fault_loss_Npm'});
writetable(T,fullfile(outDir,'pengyue_sp_tvms_healthy_crack45.csv'));

fig = figure('Visible','off','Color','w','Position',[100 100 1100 650]);
ax = axes(fig); hold(ax,'on');
plot(ax,u,kH(:,1)/1e8,'k','LineWidth',1.1);
plot(ax,u,kF(:,1)/1e8,'Color',[0.85 0.12 0.08],'LineWidth',1.0);
xlabel(ax,'Mesh phase / cycle','FontName','Times New Roman');
 ylabel(ax,'SP mesh stiffness (10^8 N/m)','FontName','Times New Roman');
legend(ax,{'Healthy','Sun-root crack, q=3 mm, gamma=45 deg'}, ...
    'Location','best','FontName','Times New Roman');
set(ax,'FontName','Times New Roman','LineWidth',0.8,'Box','on');
grid(ax,'on');
exportgraphics(fig,fullfile(outDir,'pengyue_sp_tvms_healthy_crack45.png'),'Resolution',600);
savefig(fig,fullfile(outDir,'pengyue_sp_tvms_healthy_crack45.fig'));
close(fig);

report = struct('dof',size(M,1),'meshFrequencyHz',p.kin.f_mesh, ...
    'contactRatio',modelH.contactRatio,'doublePairFraction',doubleFraction, ...
    'healthyMinNpm',min(kH(:,1)),'healthyMaxNpm',max(kH(:,1)), ...
    'healthyMeanNpm',mean(kH(:,1)),'crackQ3Gamma45MinNpm',min(kF(:,1)), ...
    'crackQ3Gamma45MaxNpm',max(kF(:,1)),'crackQ3Gamma45MeanNpm',mean(kF(:,1)), ...
    'outputCsv',fullfile(outDir,'pengyue_sp_tvms_healthy_crack45.csv'));
end
