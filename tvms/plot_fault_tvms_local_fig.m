function plot_fault_tvms_local_fig
%PLOT_FAULT_TVMS_LOCAL_FIG Create editable MATLAB .fig thesis figures.
% Plot-only operation: reads the saved q=0.50 full-recurrence CSV data and
% applies exact cyclic sample re-indexing around the maximum loss sample.
% No TVMS, fault parameter, active-pair rule, or CSV file is changed.

root = fileparts(fileparts(mfilename('fullpath')));
runDir = fullfile(root,'results','tvms','fault_tvms_20260908_local_break');
qTag = 'q050';
samplesPerCycle = 2000;
windowCycles = 6;
centerCycle = 3;

sun = readtable(fullfile(runDir,['sun_fault_sp_p1_' qTag '.csv']));
planetSP = readtable(fullfile(runDir,['planet_fault_sp_p1_' qTag '.csv']));
planetPR = readtable(fullfile(runDir,['planet_fault_pr_p1_' qTag '.csv']));

[x,sunH,sunF] = localWindow(sun,samplesPerCycle,windowCycles,centerCycle);
[xSP,spH,spF] = localWindow(planetSP,samplesPerCycle,windowCycles,centerCycle);
[xPR,prH,prF] = localWindow(planetPR,samplesPerCycle,windowCycles,centerCycle);
assert(isequal(x,xSP) && isequal(x,xPR),'Local figure phase grids differ.');

% Fig. 2-9: one representative solar-tooth fault event.
figSun = figure('Color','w','Position',[100 100 780 465], ...
    'Name','Fig2_9_sun_fault_TVMS_local','NumberTitle','off');
ax = axes(figSun); hold(ax,'on');
plot(ax,x,sunH/1e8,'Color',[0 0 0],'LineWidth',0.80);
plot(ax,x,sunF/1e8,'Color',[0.12 0.31 0.49],'LineWidth',1.35);
styleAxes(ax,'SP 啮合刚度',false);
xlabel(ax,'啮合相位（周期）','FontName','SimSun');
legend(ax,{'健康','故障'},'FontName','SimSun','FontSize',9,'Location','eastoutside','Box','off');
savefig(figSun,fullfile(runDir,'Fig2_9_sun_fault_TVMS_local.fig'));
close(figSun);

% Fig. 2-10: same planet tooth, separately centered SP and PR local views.
figPlanet = figure('Color','w','Position',[100 100 800 690], ...
    'Name','Fig2_10_planet_fault_TVMS_local','NumberTitle','off');
ax1 = subplot(2,1,1,'Parent',figPlanet); hold(ax1,'on');
plot(ax1,x,spH/1e8,'Color',[0 0 0],'LineWidth',0.80);
plot(ax1,x,spF/1e8,'Color',[0.12 0.31 0.49],'LineWidth',1.35);
styleAxes(ax1,'SP 啮合刚度',true);
legend(ax1,{'健康','故障'},'FontName','SimSun','FontSize',9,'Location','eastoutside','Box','off');

ax2 = subplot(2,1,2,'Parent',figPlanet); hold(ax2,'on');
plot(ax2,x,prH/1e8,'Color',[0 0 0],'LineWidth',0.80);
plot(ax2,x,prF/1e8,'Color',[0.65 0.29 0.00],'LineWidth',1.35);
styleAxes(ax2,'PR 啮合刚度',true);
xlabel(ax2,'啮合相位（周期）','FontName','SimSun');
legend(ax2,{'健康','故障'},'FontName','SimSun','FontSize',9,'Location','eastoutside','Box','off');

annotation(figPlanet,'textbox',[0.015 0.948 0.06 0.035],'String','(a)', ...
    'FontName','Times New Roman','FontWeight','bold','FontSize',10.5, ...
    'LineStyle','none','FitBoxToText','off');
annotation(figPlanet,'textbox',[0.015 0.482 0.06 0.035],'String','(b)', ...
    'FontName','Times New Roman','FontWeight','bold','FontSize',10.5, ...
    'LineStyle','none','FitBoxToText','off');
savefig(figPlanet,fullfile(runDir,'Fig2_10_planet_fault_TVMS_local.fig'));
close(figPlanet);

fprintf('Saved editable FIG files:\n%s\n%s\n', ...
    fullfile(runDir,'Fig2_9_sun_fault_TVMS_local.fig'), ...
    fullfile(runDir,'Fig2_10_planet_fault_TVMS_local.fig'));
end

function [x,healthy,fault] = localWindow(tbl,samplesPerCycle,windowCycles,centerCycle)
nFull = height(tbl)-1; % terminal sample duplicates the recurrence endpoint.
assert(mod(nFull,samplesPerCycle)==0,'Input CSV must have an integer number of mesh cycles.');
[maximumLoss,centerIndex] = max(tbl.delta_k_N_per_m(1:nFull));
assert(maximumLoss>0,'Input CSV contains no fault stiffness loss.');
offset = centerIndex-1-centerCycle*samplesPerCycle;
index = mod((0:windowCycles*samplesPerCycle)'+offset,nFull)+1;
x = (0:windowCycles*samplesPerCycle)'/samplesPerCycle;
healthy = tbl.k_healthy_N_per_m(index);
fault = tbl.k_fault_N_per_m(index);
assert(abs(healthy(centerCycle*samplesPerCycle+1)-tbl.k_healthy_N_per_m(centerIndex))<1e-6, ...
    'Cyclic extraction changed a healthy value.');
assert(abs(fault(centerCycle*samplesPerCycle+1)-tbl.k_fault_N_per_m(centerIndex))<1e-6, ...
    'Cyclic extraction changed a fault value.');
end

function styleAxes(ax,chineseY,showX)
set(ax,'Box','off','LineWidth',0.8,'FontName','Times New Roman','FontSize',10, ...
    'XLim',[0 6],'XTick',0:6,'TickDir','out');
xline(ax,3,'--','Color',[0.55 0.55 0.55],'LineWidth',0.65,'HandleVisibility','off');
% Chinese and unit use independent fonts, keeping all numerals/Latin unit
% characters in Times New Roman while retaining Chinese glyph coverage.
text(ax,-0.095,0.36,chineseY,'Units','normalized','Rotation',90, ...
    'FontName','SimSun','HorizontalAlignment','center','VerticalAlignment','middle','Clipping','off');
text(ax,-0.095,0.73,'(10^{8} N/m)','Units','normalized','Rotation',90, ...
    'FontName','Times New Roman','Interpreter','tex', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','Clipping','off');
if ~showX
    set(ax,'XTickLabel',[]);
end
end
