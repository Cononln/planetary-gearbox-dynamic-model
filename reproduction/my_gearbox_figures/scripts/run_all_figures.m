function run_all_figures()
%RUN_ALL_FIGURES batch runner for all 19 user-gearbox figures
    figs = {'fig3_01_sp_tvms','fig3_02_mesh_phase','fig3_03_pl_ext', ...
            'fig3_04_pl_int','fig3_05_sp31','fig3_06_pr31','fig3_07_dte_time', ...
            'fig3_08_spectra','fig3_09_depth_time','fig3_10_depth_spectrum_3d', ...
            'fig3_11_angle_time','fig3_12_angle_spectrum_3d','fig3_13_pl_time', ...
            'fig3_14_pl_spec','fig3_15_pl_zoom','fig3_16_pl_depth_time', ...
            'fig3_17_pl_depth_3d','fig3_18_pl_angle_time','fig3_19_pl_angle_3d'};
    nOK = 0;
    for k = 1:numel(figs)
        try
            feval(figs{k});
            nOK = nOK + 1;
        catch e
            disp(['FAIL ' figs{k} ': ' e.message]);
        end
    end
    disp(['completed ' num2str(nOK) '/19 figures']);
end
