function fig3_01_sp_tvms()
%Fig3_01 assembled sun-planet mesh stiffness, 2 mesh periods vs rotation
%   angle, healthy vs cracked (thesis fig 3-2 style), sun crack q1.2 g45
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(fullfile(rootDir, 'helpers'));
DATA = fullfile(rootDir, 'data');
h = readmatrix(fullfile(DATA, 'ksptime_healthy_1.csv'));
c = readmatrix(fullfile(DATA, 'ksptime_crack_1.csv'));
pitch = 2*pi/21;                          % sun tooth pitch, rad
[~, id] = min(c(:,2));                    % cracked-tooth meshing period
Tdip = c(id, 1);
T0 = max(0, min(Tdip - 1.7, 19));         % 2-period window around the dip
mh = h(:,1) >= T0 - 1e-9 & h(:,1) <= T0 + 2 + 1e-9;
mc = c(:,1) >= T0 - 1e-9 & c(:,1) <= T0 + 2 + 1e-9;
th = (h(mh,1) - T0)*pitch;
thc = (c(mc,1) - T0)*pitch;
fig = figure('Visible','off','Color','w'); ax = axes('Parent',fig);
plot(ax, th, h(mh,2), 'b-', 'LineWidth', 1.3); hold(ax,'on');
plot(ax, thc, c(mc,2), 'r--', 'LineWidth', 1.2);
xlabel(ax, '转角 \theta (rad)'); ylabel(ax, '刚度 k_{sp} (N/m)');
legend(ax, {'健康','含裂纹'}, 'Location','east');
apply_my_figure_style(fig, ax, '2d');
export_figure_all(fig, 'Fig3_01_sp_tvms', rootDir);
fprintf('done Fig3_01_sp_tvms\n');

end
