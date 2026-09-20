function fig3_07_dte_time()
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(fullfile(rootDir, 'helpers'));
DATA = fullfile(rootDir, 'data');
CB = [0 0 1]; CR = [0.85 0.10 0.10];
LAYERS = [0.47 0.47 0.47; 0.36 0.76 0.91; 0.13 0.27 0.73; 0.93 0.53 0.13; 0.80 0.13 0.13];
th = readmatrix(fullfile(DATA, 'time_healthy.csv'));
tc = readmatrix(fullfile(DATA, 'time_sun_q1.5_g45.csv'));
fig = figure('Visible','off','Color','w');
ax1 = subplot(2,1,1,'Parent',fig);
plot(ax1, th(:,1), th(:,2), 'b-', 'LineWidth', 0.5);
ylabel(ax1, 'DTE (\mum)'); ax1.Title.String = '健康'; ax1.Title.FontSize = 10;
ax2 = subplot(2,1,2,'Parent',fig);
plot(ax2, tc(:,1), tc(:,2), 'b-', 'LineWidth', 0.5);
xlabel(ax2, '时间 t (s)'); ylabel(ax2, 'DTE (\mum)');
ax2.Title.String = '太阳轮裂纹 q=1.5mm 45°'; ax2.Title.FontSize = 10;
apply_my_figure_style(fig, ax1, '2d');
export_figure_all(fig, 'Fig3_07_dte_time', rootDir);
fprintf('done Fig3_07_dte_time\n');

end
