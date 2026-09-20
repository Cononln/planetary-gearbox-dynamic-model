function fig3_17_pl_depth_3d()
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(fullfile(rootDir, 'helpers'));
DATA = fullfile(rootDir, 'data');
CB = [0 0 1]; CR = [0.85 0.10 0.10];
LAYERS = [0.47 0.47 0.47; 0.36 0.76 0.91; 0.13 0.27 0.73; 0.93 0.53 0.13; 0.80 0.13 0.13];
names = {'spec_healthy','spec_pl_q0.4_g45','spec_pl_q0.8_g45','spec_pl_q1.2_g45','spec_pl_q1.6_g45'};
L = numel(names); freq = []; amps = cell(1,L);
for k = 1:L
    d = readmatrix(fullfile(DATA, [names{k} '.csv']));
    freq = d(:,1); amps{k} = d(:,2);
end
fig = figure('Visible','off','Color','w'); ax = axes('Parent',fig);
offs = [0 0.4 0.8 1.2 1.6];
opts = struct('lineW', 0.55, 'xLim', [0 2000], 'layerGap', 1.0, 'yTickLabels', {{'0','0.4','0.8','1.2','1.6'}});
plot_my_3d_waterfall(ax, freq, amps, offs, LAYERS, opts);
xlabel(ax, '频率 f (Hz)'); ylabel(ax, '深度 q (mm)'); zlabel(ax, '幅值 (\mum)');
apply_my_figure_style(fig, ax, '3d');
export_figure_all(fig, 'Fig3_17_pl_depth_3d', rootDir);
fprintf('done Fig3_17_pl_depth_3d\n');

end
