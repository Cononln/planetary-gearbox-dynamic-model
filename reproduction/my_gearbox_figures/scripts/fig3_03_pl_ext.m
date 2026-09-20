function fig3_03_pl_ext()
%Fig3_03 sun-planet pair mesh stiffness with the cracked planet tooth as
%   the tracked meshing tooth (thesis fig 3-5 style): the crack lowers the
%   curve over the whole mesh cycle, most inside the single-tooth zone
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(fullfile(rootDir, 'helpers'));
DATA = fullfile(rootDir, 'data');
d = readmatrix(fullfile(DATA, 'tvms_pl_ext.csv'));
eta = d(:,1); kHp = d(:,2); kCp = d(:,3);   % single-pair stiffness vs eta
eps = 1.1507;                               % SP contact ratio
pitch = 2*pi/31;                            % planet tooth pitch, rad
u = linspace(0, 2, 4001);                   % two mesh periods
xi = u - floor(u);
k1H = interp1(eta, kHp, xi/eps);            % tracked (entering) tooth pair
k1C = interp1(eta, kCp, xi/eps);
m2 = xi < eps - 1;                          % double-tooth zone
k2H = interp1(eta, kHp, (1+xi)/eps);        % adjacent healthy tooth pair
k2H(~m2) = 0;                               % avoid NaN leaking via the mask
kHtot = k1H + m2.*k2H;
kCtot = k1C + m2.*k2H;
th = u*pitch;
fig = figure('Visible','off','Color','w'); ax = axes('Parent',fig);
plot(ax, th, kHtot, 'b-', 'LineWidth', 1.3); hold(ax,'on');
plot(ax, th, kCtot, 'r--', 'LineWidth', 1.2);
xlabel(ax, '转角 \theta (rad)'); ylabel(ax, '刚度 k_{sp} (N/m)');
legend(ax, {'健康','含裂纹'}, 'Location','east');
apply_my_figure_style(fig, ax, '2d');
export_figure_all(fig, 'Fig3_03_pl_ext', rootDir);
fprintf('done Fig3_03_pl_ext\n');

end
