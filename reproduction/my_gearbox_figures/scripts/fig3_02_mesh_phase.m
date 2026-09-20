function fig3_02_mesh_phase()
%Fig3_02 sun-planet mesh stiffness for the three mesh phases, 21 mesh
%   periods, healthy vs cracked (thesis fig 3-3 style), sun crack q1.2 g45
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(fullfile(rootDir, 'helpers'));
DATA = fullfile(rootDir, 'data');
fig = figure('Visible','off','Color','w');
pos3 = [0.13 0.73 0.715 0.20;       % headroom so the top ylabel is not clipped
        0.13 0.41 0.715 0.20;
        0.13 0.09 0.715 0.20];
axs = gobjects(1,3);
for i = 1:3
    hh = readmatrix(fullfile(DATA, sprintf('ksptime_healthy_%d.csv', i)));
    hc = readmatrix(fullfile(DATA, sprintf('ksptime_crack_%d.csv', i)));
    axs(i) = subplot(3,1,i,'Parent',fig);
    plot(axs(i), hh(:,1), hh(:,2), 'b-', 'LineWidth', 0.9); hold(axs(i),'on');
    plot(axs(i), hc(:,1), hc(:,2), 'r--', 'LineWidth', 0.9);
    ylabel(axs(i), sprintf('刚度 k_{sp%d} (N/m)', i));
    set(axs(i), 'XTickLabel', []);
    xlim(axs(i), [0 21]);
    axs(i).Position = pos3(i, :);
end
xlabel(axs(3), '啮合周期 t/T_m');
legend(axs(1), {'健康','含裂纹'}, 'Location','northwest');
apply_my_figure_style(fig, axs(3), '2d');
fig.Units = 'centimeters';
fig.Position = [2 2 14.52499 9.0];  % 3-row figure: taller than the 2D standard
export_figure_all(fig, 'Fig3_02_mesh_phase', rootDir);
fprintf('done Fig3_02_mesh_phase\n');

end
