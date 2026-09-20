%FIG3_10_DEPTH_SPECTRUM_3D example: 3D layered depth spectrum (user gearbox)
%   Output: Fig3_10_depth_spectrum_3d (.fig/.png/.tif/.eps/.pdf)
%   Data:   ../data/spec_healthy.csv, spec_sun_q{0.4,0.8,1.2,1.6}_g45.csv
%           (col1 f_Hz, col2 amp_um; Plan B, e_av excluded)
%   Main axes: 0-2000 Hz layered full spectrum (nf_m comb)
%   Inset:     mini-3D layered spectrum of the 0-60 Hz fault band,
%              nfs = fm/Zs = 8 Hz family (Zs=21=3x7 -> in-phase planetary,
%              single-mesh DTE fault frequency is fm/Zs, NOT Np*fm/Zs)

function fig3_10_depth_spectrum_3d()
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(fullfile(rootDir, 'helpers'));

files = {'spec_healthy.csv', 'spec_sun_q0.4_g45.csv', 'spec_sun_q0.8_g45.csv', ...
         'spec_sun_q1.2_g45.csv', 'spec_sun_q1.6_g45.csv'};
depths = [0, 0.4, 0.8, 1.2, 1.6];
L = numel(files);

freq = [];  amps = cell(1, L);
for k = 1:L
    d = readmatrix(fullfile(rootDir, 'data', files{k}));
    freq = d(:, 1);
    amps{k} = d(:, 2);
end

fig = figure('Visible', 'off', 'Color', 'w');
ax = axes('Parent', fig);

% Unified condition-to-colour mapping used by the time-domain figures:
% healthy/smallest parameter -> grey; increasing severity -> red.
colors = [0.47 0.47 0.47;
          0.36 0.76 0.91;
          0.13 0.27 0.73;
          0.93 0.53 0.13;
          0.80 0.13 0.13];

opts = struct('lineW', 0.55, 'xLim', [0 2000], ...
              'layerGap', 1.0, 'yTickLabels', ...
              {{'0', '0.4', '0.8', '1.2', '1.6'}});
plot_my_3d_waterfall(ax, freq, amps, depths, colors, opts);

xlabel(ax, '频率 f (Hz)');
ylabel(ax, '深度 q (mm)');
zlabel(ax, '幅值 (\mum)');

% ---- inset: mini-3D layered spectrum of the 0-60 Hz fault band ----
ai = axes('Parent', fig);
ai.Position = [0.065 0.60 0.30 0.34];
mskI = freq <= 60;
ampsI = cell(1, L);
for k = 1:L
    ampsI{k} = amps{k}(mskI);
end
optsI = struct('lineW', 0.55, 'xLim', [0 60], ...
               'layerGap', 1.0, 'yTickLabels', ...
               {{'0', '0.4', '0.8', '1.2', '1.6'}});
plot_my_3d_waterfall(ai, freq(mskI), ampsI, depths, colors, optsI);
xlabel(ai, '频率 f (Hz)');
ai.YLabel.String = '';
ai.ZLabel.String = '';
ai.FontSize = 7;
ai.XAxis.FontSize = 7; ai.YAxis.FontSize = 7; ai.ZAxis.FontSize = 7;

% ---- unified style + export ----
apply_my_figure_style(fig, ax, '3d');
ai.FontSize = 7; ai.XAxis.FontSize = 7; ai.YAxis.FontSize = 7; ai.ZAxis.FontSize = 7;

export_figure_all(fig, 'Fig3_10_depth_spectrum_3d', rootDir);
fprintf('done Fig3_10_depth_spectrum_3d\n');
end
