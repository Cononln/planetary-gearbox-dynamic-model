function reports = run_hybrid_validations
%RUN_HYBRID_VALIDATIONS Run the retained 18+6 elastic-ring comparison model.

setup_paths;
reports.stage1 = step1_validate_18dof;
reports.stage2 = step2_validate_tvms;
reports.stage3 = step3_validate_ring_modes;
reports.stage4 = step4_validate_paths;
reports.stage5 = step5_validate_phase_truth;
save(fullfile('results','all_hybrid_validation_reports.mat'),'reports');
fprintf('\nAll optional hybrid-model stages PASSED.\n');
end
