function export_figure_all(fig, figName, rootDir)
%EXPORT_FIGURE_ALL export .fig/.png(300dpi)/.tif(600dpi)/.eps/.pdf + log
    dirs = struct('fig', fullfile(rootDir, 'fig_matlab'), ...
                  'png', fullfile(rootDir, 'export_png'), ...
                  'tif', fullfile(rootDir, 'export_tif'), ...
                  'eps', fullfile(rootDir, 'export_eps'), ...
                  'pdf', fullfile(rootDir, 'export_pdf'));
    entries = fieldnames(dirs);
    for k = 1:numel(entries)
        if ~exist(dirs.(entries{k}), 'dir')
            mkdir(dirs.(entries{k}));
        end
    end
    fig.Renderer = 'painters';
    drawnow;
    savefig(fig, fullfile(dirs.fig, [figName '.fig']));
    for attempt = 1:3
        try
            print(fig, fullfile(dirs.png, [figName '.png']), '-dpng', '-r300');
            print(fig, fullfile(dirs.tif, [figName '.tif']), '-dtiff', '-r600');
            print(fig, fullfile(dirs.eps, [figName '.eps']), '-depsc');
            print(fig, fullfile(dirs.pdf, [figName '.pdf']), '-dpdf');
            break;
        catch ME
            if attempt == 3
                rethrow(ME);
            end
            pause(1.5);
            drawnow;
        end
    end
    fid = fopen(fullfile(rootDir, 'logs', [figName '.log']), 'w');
    fprintf(fid, 'figure: %s', figName);
    fprintf(fid, ' | exported: %s', datestr(now));
    fprintf(fid, ' | formats: fig png300 tif600 eps pdf');
    fclose(fid);
    fprintf('  exported %s (fig png tif eps pdf)', figName);
end
