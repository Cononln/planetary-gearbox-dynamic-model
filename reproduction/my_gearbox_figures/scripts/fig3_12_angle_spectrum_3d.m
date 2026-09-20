function fig3_12_angle_spectrum_3d()
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(fullfile(rootDir, 'helpers'));
DATA = fullfile(rootDir, 'data');
CB = [0 0 1]; CR = [0.85 0.10 0.10];
LAYERS = [0.47 0.47 0.47; 0.36 0.76 0.91; 0.13 0.27 0.73; 0.93 0.53 0.13; 0.80 0.13 0.13];
names = {'spec_sun_q1.2_g15','spec_sun_q1.2_g30','spec_sun_q1.2_g45','spec_sun_q1.2_g60','spec_sun_q1.2_g75'};
L = numel(names); freq = []; amps = cell(1,L);
for k = 1:L
    d = readmatrix(fullfile(DATA, [names{k} '.csv']));
    freq = d(:,1); amps{k} = d(:,2);
end
fig = figure('Visible','off','Color','w'); ax = axes('Parent',fig);
offs = [15 30 45 60 75];
opts = struct('lineW', 0.55, 'xLim', [0 2000], 'layerGap', 1.0, 'yTickLabels', {{'15','30','45','60','75'}});
plot_my_3d_waterfall(ax, freq, amps, offs, LAYERS, opts);
xlabel(ax, '频率 f (Hz)'); ylabel(ax, '角度 \gamma (°)'); zlabel(ax, '幅值 (\mum)');
apply_my_figure_style(fig, ax, '3d');
export_figure_all(fig, 'Fig3_12_angle_spectrum_3d', rootDir);
fprintf('done Fig3_12_angle_spectrum_3d\n');

end
