function fig3_09_depth_time()
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(fullfile(rootDir, 'helpers'));
DATA = fullfile(rootDir, 'data');
CB = [0 0 1]; CR = [0.85 0.10 0.10];
LAYERS = [0.47 0.47 0.47; 0.36 0.76 0.91; 0.13 0.27 0.73; 0.93 0.53 0.13; 0.80 0.13 0.13];
names = {'time_healthy','time_sun_q0.4_g45','time_sun_q0.8_g45','time_sun_q1.2_g45','time_sun_q1.6_g45'};
L = numel(names); t = []; amps = cell(1,L);
for k = 1:L
    d = readmatrix(fullfile(DATA, [names{k} '.csv']));
    m = d(:,1) <= 1; t = d(m,1); amps{k} = d(m,2);   % 1 s window
    amps{k} = amps{k} - mean(amps{k});
end
fig = figure('Visible','off','Color','w'); ax = axes('Parent',fig);
offs = [0 0.4 0.8 1.2 1.6];
opts = struct('lineW', 0.40, 'mode', 'time', 'xLim', [0 1], ...
    'view', [-37.5 30], 'timeTicks', 0:0.25:1, 'timePointLimit', 16000, ...
    'preserveExtrema', true, ...
    'yTickLabels', {{'0','0.4','0.8','1.2','1.6'}});
plot_my_3d_waterfall(ax, t, amps, offs, LAYERS, opts);
xlabel(ax, '深度 q (mm)'); ylabel(ax, '时间 t (s)');
zlabel(ax, '系统响应 (\mum)');
apply_my_figure_style(fig, ax, '3d');
export_figure_all(fig, 'Fig3_09_depth_time_raw_paired', rootDir);
fprintf('done Fig3_09_depth_time_raw_paired\n');

end
