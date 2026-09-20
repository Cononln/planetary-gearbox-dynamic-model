function fig3_14_pl_spec()
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(fullfile(rootDir, 'helpers'));
DATA = fullfile(rootDir, 'data');
CB = [0 0 1]; CR = [0.85 0.10 0.10];
LAYERS = [0.47 0.47 0.47; 0.36 0.76 0.91; 0.13 0.27 0.73; 0.93 0.53 0.13; 0.80 0.13 0.13];
d = readmatrix(fullfile(DATA, 'spec_pl_q1.5_g45.csv'));
fig = figure('Visible','off','Color','w'); ax = axes('Parent',fig);
plot(ax, d(:,1), d(:,2), 'b-', 'LineWidth', 0.6);
xlim(ax, [0 2000]); xlabel(ax, '频率 f (Hz)'); ylabel(ax, '幅值 (\mum)');
apply_my_figure_style(fig, ax, '2d');
export_figure_all(fig, 'Fig3_14_pl_spec', rootDir);
fprintf('done Fig3_14_pl_spec\n');

end
