function fig3_15_pl_zoom()
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(fullfile(rootDir, 'helpers'));
DATA = fullfile(rootDir, 'data');
CB = [0 0 1]; CR = [0.85 0.10 0.10];
LAYERS = [0.47 0.47 0.47; 0.36 0.76 0.91; 0.13 0.27 0.73; 0.93 0.53 0.13; 0.80 0.13 0.13];
d = readmatrix(fullfile(DATA, 'spec_pl_q1.5_g45.csv'));
fp = 168/31; fm = 168;
fig = figure('Visible','off','Color','w');
ax1 = subplot(2,1,1,'Parent',fig);
m1 = d(:,1) <= 120;
plot(ax1, d(m1,1), d(m1,2), 'b-', 'LineWidth', 0.8); hold(ax1,'on');
for k = 1:10
    f0 = k*fp; m2 = abs(d(:,1)-f0) <= 1.0;
    plot(ax1, d(m2,1), d(m2,2), 'r.', 'MarkerSize', 7);
end
xlabel(ax1, '频率 f (Hz)'); ylabel(ax1, '幅值 (\mum)');
ax1.Title.String = '低频放大 (n_f_p = 5.42 Hz)'; ax1.Title.FontSize = 10;
ax2 = subplot(2,1,2,'Parent',fig);
m3 = abs(d(:,1)-6*fm) <= 80;
plot(ax2, d(m3,1), d(m3,2), 'b-', 'LineWidth', 0.8); hold(ax2,'on');
for sg = [-1 1]
    for k = 1:3
        f0 = 6*fm + sg*k*fp; m2 = abs(d(:,1)-f0) <= 1.0;
        plot(ax2, d(m2,1), d(m2,2), 'g.', 'MarkerSize', 7);
    end
end
xlim(ax2, [6*fm-80 6*fm+80]); ylim(ax2, [0 0.001]);
xlabel(ax2, '频率 f (Hz)'); ylabel(ax2, '幅值 (\mum)');
ax2.Title.String = '6f_m \pm n_f_p 放大'; ax2.Title.FontSize = 10;
apply_my_figure_style(fig, ax1, '2d');
export_figure_all(fig, 'Fig3_15_pl_zoom', rootDir);
fprintf('done Fig3_15_pl_zoom\n');

end
