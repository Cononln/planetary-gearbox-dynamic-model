function result = diagnose_fe_modal_path_gain(options)
%DIAGNOSE_FE_MODAL_PATH_GAIN Export contact-angle/modal path attenuation.
%
% Input: mass-normalized ANSYS ring asset and frequency/modal-band settings.
% Output: complex contact-to-sensor transfer and modal-energy diagnostics
% versus moving PR contact angle, saved under OutputDir. Source dynamics
% and prior result directories are not changed.

arguments
    options.AssetPath (1,1) string = ""
    options.OutputDir (1,1) string = ""
    options.ModalBandHz (1,2) double {mustBeNonnegative} = [15000 30000]
    options.ModalDampingRatio (1,1) double {mustBePositive} = 0.02
    options.NAngle (1,1) double {mustBeInteger,mustBeGreaterThan(options.NAngle,7)} = 720
    options.FrequencyHz (:,1) double = [168; 24; 336]
end
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_paths;
if strlength(options.AssetPath)==0
    assetPath = fullfile(projectRoot,'results','fe_ring_asset', ...
        'ansys_ring_modal_transfer_asset_exact_four_sensor.mat');
else
    assetPath = char(options.AssetPath);
end
model = pg_fe_ring_modal_model(string(assetPath), ...
    'ModalDampingRatio',options.ModalDampingRatio);
keep = model.frequencyHz >= options.ModalBandHz(1) & ...
    model.frequencyHz <= options.ModalBandHz(2);
if ~any(keep)
    error('diagnose_fe_modal_path_gain:NoModes', ...
        'No modes fall in the requested modal band.');
end
model.frequencyHz = model.frequencyHz(keep);
model.omega = model.omega(keep);
model.zeta = model.zeta(keep);
model.inputShapeNormal = model.inputShapeNormal(:,keep);
model.sensorShapeRadial = model.sensorShapeRadial(:,keep);
model.retainedModeIndex = find(keep);

if strlength(options.OutputDir)==0
    stamp = datestr(now,'yyyymmdd_HHMMSS');
    outDir = fullfile(projectRoot,'results', ...
        ['fe_modal_path_gain_diagnostic_' stamp]);
else
    outDir = char(options.OutputDir);
end
if ~exist(outDir,'dir'); mkdir(outDir); end

theta = (0:options.NAngle-1)'*(2*pi/options.NAngle);
[H,E,detail] = pg_fe_ring_modal_path_gain(theta,options.FrequencyHz(:),model);
sensorAngles = model.sensorAnglesLocalDeg(:).';

% Flatten a machine-readable complex/energy table.
nT = numel(theta); nS = numel(sensorAngles); nF = numel(options.FrequencyHz);
rows = nT*nS*nF;
thetaCol = zeros(rows,1); sensorCol=zeros(rows,1); freqCol=zeros(rows,1);
realCol=zeros(rows,1); imagCol=zeros(rows,1); magCol=zeros(rows,1);
phaseCol=zeros(rows,1); energyCol=zeros(rows,1);
r=0;
for jf=1:nF
    for is=1:nS
        ii=r+(1:nT); r=r+nT;
        thetaCol(ii)=theta; sensorCol(ii)=sensorAngles(is); freqCol(ii)=options.FrequencyHz(jf);
        z=H(:,is,jf); realCol(ii)=real(z); imagCol(ii)=imag(z);
        magCol(ii)=abs(z); phaseCol(ii)=rad2deg(angle(z));
        energyCol(ii)=E(:,is,jf);
    end
end
T=table(thetaCol,rad2deg(thetaCol),sensorCol,freqCol,realCol,imagCol, ...
    magCol,phaseCol,energyCol, ...
    'VariableNames',{'contact_angle_rad','contact_angle_deg','sensor_angle_deg', ...
    'frequency_Hz','H_real_per_N','H_imag_per_N','H_abs_per_N', ...
    'H_phase_deg','energy_gain_per_N2'});
writetable(T,fullfile(outDir,'contact_sensor_modal_path_gain.csv'));

% Angle-wise relative gain and integrated band energy (no distance model).
bandF = linspace(options.ModalBandHz(1),options.ModalBandHz(2),121).';
[Hb,Eb] = pg_fe_ring_modal_path_gain(theta,bandF,model);
bandEnergy = squeeze(trapz(bandF,Eb,3));
bandMeanSquare = bandEnergy/(bandF(end)-bandF(1));
bandRmsGain = sqrt(max(bandMeanSquare,0));
Tb=table(theta,rad2deg(theta),'VariableNames',{'contact_angle_rad','contact_angle_deg'});
for is=1:nS
    token=sprintf('sensor_%03ddeg',round(sensorAngles(is)));
    Tb.([token '_band_energy'])=bandEnergy(:,is);
    Tb.([token '_band_mean_square_gain'])=bandMeanSquare(:,is);
    Tb.([token '_band_rms_gain'])=bandRmsGain(:,is);
end
writetable(Tb,fullfile(outDir,'contact_sensor_modal_band_energy.csv'));

% Long table adds source/sensor angular separation for diagnosis only.  The
% separation is never used to calculate the modal transfer or energy.
nLong=nT*nS;
contactLong=zeros(nLong,1); localLong=zeros(nLong,1); globalLong=zeros(nLong,1);
distanceLong=zeros(nLong,1); energyLong=zeros(nLong,1); rmsLong=zeros(nLong,1);
normalizedLong=zeros(nLong,1); r=0;
for is=1:nS
    ii=r+(1:nT); r=r+nT;
    contactLong(ii)=rad2deg(theta);
    localLong(ii)=sensorAngles(is);
    globalLong(ii)=model.sensorAnglesGlobalDeg(is);
    distanceLong(ii)=abs(mod(contactLong(ii)-globalLong(ii)+180,360)-180);
    energyLong(ii)=bandEnergy(:,is);
    rmsLong(ii)=bandRmsGain(:,is);
    normalizedLong(ii)=bandEnergy(:,is)/max(bandEnergy(:,is));
end
Tlong=table(contactLong,localLong,globalLong,distanceLong,energyLong,rmsLong,normalizedLong, ...
    'VariableNames',{'contact_angle_deg','sensor_local_angle_deg', ...
    'sensor_global_angle_deg','angular_distance_deg','band_energy', ...
    'band_rms_gain','energy_normalized_by_sensor_max'});
writetable(Tlong,fullfile(outDir,'contact_sensor_modal_band_energy_long.csv'));

% Constraint/path asymmetry summary.  Near/far regions are descriptive
% bins only, not a propagation law or correction factor.
meanEnergy=zeros(nS,1); maxEnergy=zeros(nS,1); maxAngle=zeros(nS,1);
nearMean=zeros(nS,1); farMean=zeros(nS,1); nearFarRatio=zeros(nS,1);
distanceCorr=zeros(nS,1);
for is=1:nS
    distance=abs(mod(rad2deg(theta)-model.sensorAnglesGlobalDeg(is)+180,360)-180);
    energy=bandEnergy(:,is);
    meanEnergy(is)=mean(energy); [maxEnergy(is),ix]=max(energy);
    maxAngle(is)=rad2deg(theta(ix));
    nearMean(is)=mean(energy(distance<=15));
    farMean(is)=mean(energy(distance>=150));
    nearFarRatio(is)=nearMean(is)/farMean(is);
    R=corrcoef(distance,energy); distanceCorr(is)=R(1,2);
end
Ts=table(sensorAngles(:),model.sensorAnglesGlobalDeg(:),meanEnergy,maxEnergy, ...
    maxAngle,nearMean,farMean,nearFarRatio,distanceCorr, ...
    'VariableNames',{'sensor_local_angle_deg','sensor_global_angle_deg', ...
    'mean_band_energy','max_band_energy','max_energy_contact_angle_deg', ...
    'near_15deg_mean_energy','far_150deg_mean_energy', ...
    'near_to_far_energy_ratio','corr_angular_distance_energy'});
writetable(Ts,fullfile(outDir,'sensor_modal_path_energy_summary.csv'));

% Per-mode participation/energy summaries, useful for constraints diagnosis.
modeEnergy=zeros(numel(model.frequencyHz),nS);
for im=1:numel(model.frequencyHz)
    % Integrate source shape energy over contact angle and retain sensor shape.
    sourceRms=sqrt(mean(model.inputShapeNormal(:,im).^2));
    modeEnergy(im,:) = sourceRms^2 * model.sensorShapeRadial(:,im).'.^2;
end
Tm=table(model.retainedModeIndex(:),model.frequencyHz(:), ...
    'VariableNames',{'full_mode_index','mode_frequency_Hz'});
for is=1:nS
    Tm.(sprintf('sensor_%03ddeg_mode_energy',round(sensorAngles(is))))=modeEnergy(:,is);
end
writetable(Tm,fullfile(outDir,'modal_energy_by_sensor.csv'));

% Minimal diagnostic plot: no response data or artificial scaling.
fig=figure('Color','w','Visible','off','Position',[100 100 1200 760]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for jf=1:nF
    nexttile; hold on;
    for is=1:nS
        plot(rad2deg(theta),20*log10(max(abs(H(:,is,jf)),realmin)), ...
            'LineWidth',1.0,'DisplayName',sprintf('%g deg',sensorAngles(is)));
    end
    grid on; box on; xlabel('Contact angle / deg'); ylabel('|H| / dB re 1 N');
    title(sprintf('Modal path magnitude at %.6g Hz',options.FrequencyHz(jf)));
    if jf==1; legend('Location','best'); end
end
nexttile; hold on;
for is=1:nS
    plot(rad2deg(theta),10*log10(max(bandEnergy(:,is),realmin)), ...
        'LineWidth',1.0,'DisplayName',sprintf('%g deg',sensorAngles(is)));
end
grid on; box on; xlabel('Contact angle / deg'); ylabel('Integrated modal energy / dB');
title(sprintf('Integrated modal band %.0f–%.0f Hz',options.ModalBandHz));
legend('Location','best');
exportgraphics(fig,fullfile(outDir,'Fig_modal_path_gain_vs_contact_angle.png'),'Resolution',300);
exportgraphics(fig,fullfile(outDir,'Fig_modal_path_gain_vs_contact_angle.pdf'),'ContentType','vector');
close(fig);

fid=fopen(fullfile(outDir,'README.txt'),'w');
fprintf(fid,'ANSYS modal-shape path-gain diagnostic\n');
fprintf(fid,'Asset: %s\n',assetPath);
fprintf(fid,'Retained modes: %d/%d, frequency range %.9g..%.9g Hz\n',sum(keep),numel(keep),min(model.frequencyHz),max(model.frequencyHz));
fprintf(fid,'Sensors local angles (deg): %s\n',mat2str(sensorAngles));
fprintf(fid,'Damping ratio: %.9g (provisional; not exported by ANSYS)\n',options.ModalDampingRatio);
fprintf(fid,'H_s(theta,w)=sum_m psi_in(theta,m)*Hacc_m(w)*psi_sensor(s,m).\n');
fprintf(fid,'Energy gain is |H|^2 per N^2; no empirical distance weighting or fitted gains.\n');
fprintf(fid,'This diagnostic does not claim a causal arrival time. Use a time-domain modal state solver for that.\n');
fclose(fid);

result=struct();
result.outputDir=outDir;
result.model=model;
result.thetaRad=theta;
result.frequencyHz=options.FrequencyHz(:);
result.H=H;
result.energy=E;
result.bandFrequencyHz=bandF;
result.bandEnergy=bandEnergy;
result.tables={T,Tb,Tlong,Ts,Tm};
fprintf('Modal path-gain diagnostic written to %s\n',outDir);
end
