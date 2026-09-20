function fig3_08_spectra()
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(fullfile(rootDir, 'helpers'));
DATA = fullfile(rootDir, 'data');
CB = [0 0 1]; CR = [0.85 0.10 0.10];
LAYERS = [0.47 0.47 0.47; 0.36 0.76 0.91; 0.13 0.27 0.73; 0.93 0.53 0.13; 0.80 0.13 0.13];
sh = readmatrix(fullfile(DATA, 'spec_healthy.csv'));
sc = readmatrix(fullfile(DATA, 'spec_sun_q1.5_g45.csv'));
fm = 168; fs = fm/21;
fig = figure('Visible','off','Color','w');
ax1 = subplot(2,2,1,'Parent',fig);
plot(ax1, sh(:,1), sh(:,2), 'b-', 'LineWidth', 0.6);
xlim(ax1, [0 2000]); xlabel(ax1,'频率 f (Hz)'); ylabel(ax1,'幅值 (\mum)');
ax1.Title.String = '健康'; ax1.Title.FontSize = 10;
ax2 = subplot(2,2,2,'Parent',fig);
plot(ax2, sc(:,1), sc(:,2), 'b-', 'LineWidth', 0.6);
xlim(ax2, [0 2000]); xlabel(ax2,'频率 f (Hz)'); ylabel(ax2,'幅值 (\mum)');
ax2.Title.String = '太阳轮裂纹'; ax2.Title.FontSize = 10;
ax3 = subplot(2,2,3,'Parent',fig);
m3 = sc(:,1) <= 160;
plot(ax3, sc(m3,1), sc(m3,2), 'b-', 'LineWidth', 0.8); hold(ax3,'on');
for k = 1:6
    f0 = k*fs; m2 = abs(sc(:,1)-f0) <= 1.0;
    [am, ii] = max(sc(m2,2));
    plot(ax3, sc(m2,1), sc(m2,2), 'r.', 'MarkerSize', 7);
end
xlim(ax3, [0 160]); xlabel(ax3,'频率 f (Hz)'); ylabel(ax3,'幅值 (\mum)');
ax3.Title.String = '低频放大 (n_f_s = 8 Hz)'; ax3.Title.FontSize = 10;
ax4 = subplot(2,2,4,'Parent',fig);
m4 = abs(sc(:,1)-6*fm) <= 80;
plot(ax4, sc(m4,1), sc(m4,2), 'b-', 'LineWidth', 0.8); hold(ax4,'on');
for sg = [-1 1]
    for k = 1:3
        f0 = 6*fm + sg*k*fs; m2 = abs(sc(:,1)-f0) <= 1.0;
        plot(ax4, sc(m2,1), sc(m2,2), 'g.', 'MarkerSize', 7);
    end
end
xlim(ax4, [6*fm-80 6*fm+80]); ylim(ax4, [0 0.001]);
xlabel(ax4,'频率 f (Hz)'); ylabel(ax4,'幅值 (\mum)');
ax4.Title.String = '6f_m \pm n_f_s 放大'; ax4.Title.FontSize = 10;
apply_my_figure_style(fig, ax1, '2d');
export_figure_all(fig, 'Fig3_08_spectra', rootDir);
fprintf('done Fig3_08_spectra\n');

end
