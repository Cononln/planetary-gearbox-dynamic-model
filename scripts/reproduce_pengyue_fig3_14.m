function output = reproduce_pengyue_fig3_14(options)
%REPRODUCE_PENGYUE_FIG3_14 Reproduce the five-angle 3-D spectrum in Fig. 3-14.
% The physical result is saved in micrometre DTE amplitude. A second figure
% uses one GLOBAL display scale so the z-axis matches the original paper's
% 0-0.02 visual range. That global scale is a plotting aid only; it is not
% evidence that the source paper's undisclosed displacement normalization bc
% has been recovered.

arguments
    options.Visible (1,1) string {mustBeMember(options.Visible,["on","off"])} = "on"
    options.RecordCarrierCycles (1,1) double {mustBeInteger,mustBePositive} = 17
    options.DiscardCarrierCycles (1,1) double {mustBeInteger,mustBeNonnegative} = 2
    options.RelTol (1,1) double {mustBePositive} = 1e-7
    options.AbsTol (1,1) double {mustBePositive} = 1e-9
    options.PaperStylePeak (1,1) double {mustBePositive} = 0.018
end

setup_paths;
validation = validate_pengyue_fig3_14;
p = pg_parameters_pengyue_fig314;
angles = p.reproduction.fig314.crackAnglesDeg;
qMm = p.reproduction.fig314.crackDepthMm;

carrierPeriod = 1/p.kin.f_c;
discardTime = options.DiscardCarrierCycles*carrierPeriod;
recordDuration = options.RecordCarrierCycles*carrierPeriod;
totalDuration = discardTime+recordDuration;

nCase = numel(angles);
responses = cell(1,nCase);
for ic = 1:nCase
    fprintf('Simulating Peng Yue Fig. 3-14: q=%.1f mm, gamma=%g deg ...\n', ...
        qMm,angles(ic));
    fault = pg_pengyue_root_crack_case(qMm,angles(ic),p);
    r = pg_simulate_pengyue_dte(totalDuration,p,fault, ...
        'fs',p.model.fs,'RelTol',options.RelTol,'AbsTol',options.AbsTol);
    keep = r.time>=discardTime;
    r.time = r.time(keep)-discardTime;
    r.q = r.q(keep,:);
    r.qd = r.qd(keep,:);
    r.dteSunPlanet = r.dteSunPlanet(keep,:);
    r.transmissionErrorSunPlanet = r.transmissionErrorSunPlanet(keep,:);
    responses{ic} = r;
end

% Use x_sp^1, matching the thesis definition of the system DTE response.
time = responses{1}.time;
n = numel(time);
frequency = (0:floor(n/2))'*p.model.fs/n;
spectrumM = zeros(numel(frequency),nCase);
for ic = 1:nCase
    x = responses{ic}.dteSunPlanet(:,1);
    x = x-mean(x);
    X = fft(x);
    a = 2*abs(X(1:numel(frequency)))/n;
    a(1) = abs(X(1))/n;
    spectrumM(:,ic) = a;
end
spectrumUm = 1e6*spectrumM;

inMain = frequency<=p.reproduction.fig314.maxFrequencyHz;
physicalPeak = max(spectrumM(inMain,:),[],'all');
if physicalPeak<=0 || ~isfinite(physicalPeak)
    error('reproduce_pengyue_fig3_14:ZeroSpectrum','Invalid DTE spectrum.');
end
globalPaperScale = options.PaperStylePeak/physicalPeak;
spectrumPaper = spectrumM*globalPaperScale;

colors = [0.36 0.36 0.36; ...
          0.43 0.78 0.94; ...
          0.05 0.28 0.78; ...
          0.95 0.47 0.08; ...
          0.90 0.12 0.12];

if ~exist('results','dir'); mkdir('results'); end
physicalFile = fullfile('results','pengyue_fig3_14_physical.png');
paperStyleFile = fullfile('results','pengyue_fig3_14_paper_style.png');

makeFigure(spectrumUm,'DTE amplitude (\mum)',[],physicalFile,false);
makeFigure(spectrumPaper,'幅值',[0 0.02],paperStyleFile,true);

% Long-format export for independent checking and replotting.
maskExport = frequency<=p.reproduction.fig314.maxFrequencyHz;
rows = nnz(maskExport)*nCase;
angleColumn = zeros(rows,1);
frequencyColumn = zeros(rows,1);
physicalColumn = zeros(rows,1);
paperColumn = zeros(rows,1);
row = 0;
for ic = 1:nCase
    idx = find(maskExport);
    rr = row+(1:numel(idx));
    angleColumn(rr) = angles(ic);
    frequencyColumn(rr) = frequency(idx);
    physicalColumn(rr) = spectrumUm(idx,ic);
    paperColumn(rr) = spectrumPaper(idx,ic);
    row = row+numel(idx);
end
spectraTable = table(angleColumn,frequencyColumn,physicalColumn,paperColumn, ...
    'VariableNames',{'crack_angle_deg','frequency_Hz', ...
    'dte_amplitude_um','paper_style_global_scaled_amplitude'});
writetable(spectraTable,fullfile('results','pengyue_fig3_14_spectra.csv'));

summary = table(angles.',repmat(qMm,nCase,1), ...
    repmat(p.kin.f_mesh,nCase,1),repmat(p.kin.f_sun_fault,nCase,1), ...
    max(spectrumUm(inMain,:),[],1).', ...
    'VariableNames',{'crack_angle_deg','crack_depth_mm','mesh_frequency_Hz', ...
    'sun_fault_frequency_Hz','maximum_DTE_spectrum_um'});
writetable(summary,fullfile('results','pengyue_fig3_14_summary.csv'));

output.parameters = p;
output.validation = validation;
output.anglesDeg = angles;
output.crackDepthMm = qMm;
output.responses = responses;
output.frequencyHz = frequency;
output.spectrumDteM = spectrumM;
output.spectrumDteUm = spectrumUm;
output.spectrumPaperStyle = spectrumPaper;
output.paperStyleGlobalScale = globalPaperScale;
output.recordDuration = recordDuration;
output.discardTime = discardTime;
output.summary = summary;
save(fullfile('results','pengyue_fig3_14_reproduction.mat'),'output','-v7.3');

fprintf('\nPeng Yue Fig. 3-14 reproduction finished.\n');
fprintf('  physical figure   : %s\n',physicalFile);
fprintf('  paper-style figure: %s\n',paperStyleFile);
fprintf('  fm / fs_fault     : %.3f / %.3f Hz\n',p.kin.f_mesh,p.kin.f_sun_fault);
fprintf('  record duration   : %.6f s (%d carrier cycles)\n', ...
    recordDuration,options.RecordCarrierCycles);
fprintf('  paper display scale: %.6g (single global factor)\n',globalPaperScale);
fprintf('  NOTE: paper-style amplitude scaling is visual only.\n');

    function makeFigure(amplitude,zLabel,zLimits,fileName,isPaperStyle)
        fig = figure('Visible',char(options.Visible),'Color','w', ...
            'Units','pixels','Position',[80 80 980 610]);
        ax = axes(fig,'Position',[0.09 0.12 0.84 0.80]); hold(ax,'on');
        for jc = 1:nCase
            plot3(ax,frequency(inMain),angles(jc)*ones(nnz(inMain),1), ...
                amplitude(inMain,jc),'Color',colors(jc,:),'LineWidth',1.05);
        end
        xlim(ax,[0 p.reproduction.fig314.maxFrequencyHz]);
        ylim(ax,[15 75]); yticks(ax,angles);
        if ~isempty(zLimits); zlim(ax,zLimits); end
        xlabel(ax,'频率 f (Hz)','FontName','Microsoft YaHei','FontSize',11);
        ylabel(ax,'角度 \gamma (°)','FontName','Microsoft YaHei','FontSize',11);
        zlabel(ax,zLabel,'FontName','Microsoft YaHei','FontSize',11);
        set(ax,'FontName','Times New Roman','FontSize',10,'LineWidth',0.8, ...
            'Box','on','GridAlpha',0.16,'Projection','perspective');
        grid(ax,'on'); view(ax,[-38 24]);

        % Low-frequency inset highlighting n*f_s, as in the thesis figure.
        ax2 = axes(fig,'Position',[0.18 0.61 0.31 0.26]); hold(ax2,'on');
        low = frequency<=p.reproduction.fig314.lowFrequencyMaxHz;
        for jc = 1:nCase
            plot3(ax2,frequency(low),angles(jc)*ones(nnz(low),1), ...
                amplitude(low,jc),'Color',colors(jc,:),'LineWidth',0.85);
        end
        xlim(ax2,[0 p.reproduction.fig314.lowFrequencyMaxHz]);
        ylim(ax2,[15 75]); yticks(ax2,angles);
        set(ax2,'FontName','Times New Roman','FontSize',8,'LineWidth',0.65, ...
            'Box','on','GridAlpha',0.12,'Projection','perspective');
        grid(ax2,'on'); view(ax2,[-38 24]);
        if isPaperStyle
            zlim(ax2,[0 max(0.0035,0.22*options.PaperStylePeak)]);
        end
        fMark = 3*p.kin.f_sun_fault;
        [~,ii] = min(abs(frequency-fMark));
        zMark = max(amplitude(ii,:));
        text(ax2,fMark,72,1.10*zMark,'n f_s','Interpreter','tex', ...
            'FontName','Times New Roman','FontSize',10,'FontAngle','italic');

        exportgraphics(fig,fileName,'Resolution',300);
        savefig(fig,strrep(fileName,'.png','.fig'));
        if options.Visible=="off"; close(fig); end
    end
end
