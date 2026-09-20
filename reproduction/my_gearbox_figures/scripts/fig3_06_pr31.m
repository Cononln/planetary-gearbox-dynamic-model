function fig3_06_pr31()
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(fullfile(rootDir, 'helpers'));
DATA = fullfile(rootDir, 'data');
CB = [0 0 1]; CR = [0.85 0.10 0.10];
LAYERS = [0.47 0.47 0.47; 0.36 0.76 0.91; 0.13 0.27 0.73; 0.93 0.53 0.13; 0.80 0.13 0.13];
h = readmatrix(fullfile(DATA, 'kpr31_healthy.csv'));
c = readmatrix(fullfile(DATA, 'kpr31_crack.csv'));
fig = figure('Visible','off','Color','w'); ax = axes('Parent',fig);
plot(ax, h(:,1), h(:,2), 'b-', 'LineWidth', 0.8); hold(ax,'on');
plot(ax, c(:,1), c(:,2), 'r--', 'LineWidth', 0.8);
xlim(ax, [0 31]);
xlabel(ax, '啮合周期 t/T_m'); ylabel(ax, '刚度 k_{pr} (N/m)');
legend(ax, {'健康','含裂纹'}, 'Location','southeast');
apply_my_figure_style(fig, ax, '2d');
export_figure_all(fig, 'Fig3_06_pr31', rootDir);
fprintf('done Fig3_06_pr31\n');

end
