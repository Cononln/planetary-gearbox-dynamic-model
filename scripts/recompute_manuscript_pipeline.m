function status = recompute_manuscript_pipeline(includeSeverity)
%RECOMPUTE_MANUSCRIPT_PIPELINE Recompute all numerical results used by the paper.
% MATLAB is used only for simulation, estimation, diagnostics, and source
% data export. Publication figures are rendered separately in Python.

if nargin < 1
    includeSeverity = true;
end

setup_paths;
rootDir = fileparts(fileparts(mfilename('fullpath')));
originalDir = pwd;
cleanup = onCleanup(@() cd(originalDir)); %#ok<NASGU>
cd(rootDir);

if ~exist('results','dir')
    mkdir('results');
end

stageName = strings(0,1);
startedAt = NaT(0,1);
elapsedSec = zeros(0,1);
passed = false(0,1);
message = strings(0,1);

    function runStage(name,fun)
        fprintf('\n===== %s =====\n',name);
        t0 = datetime('now');
        timer = tic;
        ok = false;
        note = "";
        try
            fun();
            ok = true;
        catch exception
            note = string(getReport(exception,'extended','hyperlinks','off'));
        end
        stageName(end+1,1) = string(name);
        startedAt(end+1,1) = t0;
        elapsedSec(end+1,1) = toc(timer);
        passed(end+1,1) = ok;
        message(end+1,1) = note;
        status = table(stageName,startedAt,elapsedSec,passed,message);
        writetable(status,fullfile('results', ...
            'manuscript_recompute_status.csv'));
        if ~ok
            error('recompute_manuscript_pipeline:StageFailed', ...
                'Stage "%s" failed:\n%s',name,note);
        end
    end

runStage('Hybrid dynamic-model validation',@() run_hybrid_validations());
runStage('Root-crack TVMS validation',@() validate_root_crack_model());
runStage('Corrected 20-Nm root-crack responses', ...
    @() simulate_root_crack_responses());
runStage('Moving-path amplitude/phase truth', ...
    @() analyze_loaded_path_am_pm());
runStage('Synthetic tacholess joint phase demo', ...
    @() run_tacholess_joint_phase_demo());
runStage('Measured dual-channel tacholess phase', ...
    @() run_real_tacholess_joint_phase());
runStage('Measured phase-aligned TPSVD diagnosis', ...
    @() run_real_phase_aligned_diagnosis());

healthyFile = fullfile(rootDir,'data_local','bl','healthy_record.MAT');
pf50File = fullfile(rootDir,'data_local','planet','pf50_record.MAT');
runStage('Three-channel healthy phase decomposition', ...
    @() run_real_three_channel_phase(healthyFile));
runStage('Three-channel PF50 phase decomposition', ...
    @() run_real_three_channel_phase(pf50File));
runStage('Three-channel safe TPSVD diagnosis', ...
    @() run_three_channel_safe_diagnosis());

if includeSeverity
    severityFiles = {
        fullfile(rootDir,'data_local','planet','P025S600L00.tdms.mat')
        fullfile(rootDir,'data_local','planet','P050S600L00.tdms.mat')
        fullfile(rootDir,'data_local','planet','P075S600L00.tdms.mat')
        fullfile(rootDir,'data_local','sun','Sun050S600L00.tdms.mat')};
    severityNames = {
        'Three-channel P025 severity phase decomposition'
        'Three-channel P050 severity phase decomposition'
        'Three-channel P075 severity phase decomposition'
        'Three-channel Sun050 phase decomposition'};
    for k = 1:numel(severityFiles)
        file = severityFiles{k};
        runStage(severityNames{k},@() run_real_three_channel_phase(file));
    end
end

runStage('Controlled six-configuration phase ablation', ...
    @() run_real_phase_ablation(healthyFile,pf50File));
runStage('HMPC/PCC method comparison', ...
    @() run_hmpc_pcc_comparison(healthyFile,pf50File));
runStage('MSSP figure source-data export', ...
    @() export_mssp_core_figure_data());
runStage('MATLAB unit and regression tests',@runTests);

status = table(stageName,startedAt,elapsedSec,passed,message);
writetable(status,fullfile('results','manuscript_recompute_status.csv'));
save(fullfile('results','manuscript_recompute_status.mat'),'status');
fprintf('\nAll manuscript numerical stages PASSED.\n');
end

function runTests
projectRoot = fileparts(fileparts(mfilename('fullpath')));
suite = testsuite(fullfile(projectRoot,'tests'), ...
    'IncludeSubfolders',true);
result = run(suite);
disp(table(result));
assert(all([result.Passed]),'At least one MATLAB test failed.');
end
