function reports = run_all_validations
%RUN_ALL_VALIDATIONS Run the default lumped-parameter validations.

setup_paths;
reports.stage1 = step1_validate_18dof;
reports.stage2 = step2_validate_tvms;
reports.stage3 = step3_validate_lpm_response;
save(fullfile('results','all_lpm_validation_reports.mat'),'reports');
fprintf('\nAll default LPM stages PASSED.\n');
end
