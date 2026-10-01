function result = run_causal_fe_ring_path_comparison(options)
%RUN_CAUSAL_FE_RING_PATH_COMPARISON Compare healthy/cracked FE-path response.
% Input: a mass-normalized ANSYS ring asset and one 18-DOF source model.
% Output: signed sensor time series, event metrics, diagnostic figures and
% a README under OutputDir. All acceleration values remain in m/s^2 in CSV.
% This entry point does not add FE coordinates to the 18-DOF source system.

arguments
    options.AssetPath (1,1) string = ""
    options.ModalDampingRatio (1,1) double {mustBePositive} = 0.02
    options.RecordDuration (1,1) double {mustBePositive} = 2.0
    options.DiscardDuration (1,1) double {mustBeNonnegative} = 0.2
    options.fs (1,1) double {mustBePositive} = 524288
    options.CrackDepthMm (1,1) double {mustBeNonnegative} = 0.90
    options.CrackAngleDeg (1,1) double = 45
    options.OutputDir (1,1) string = ""
end
%% 1. Resolve paths and model inputs
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_paths;
if strlength(options.AssetPath)==0
    assetDir = fullfile(projectRoot,'results','fe_ring_asset');
    fourSensorPath = fullfile(assetDir, ...
        'ansys_ring_modal_transfer_asset_exact_four_sensor.mat');
    threeSensorPath = fullfile(assetDir, ...
        'ansys_ring_modal_transfer_asset_exact_three_sensor.mat');
    if isfile(fourSensorPath)
        assetPath = fourSensorPath;
    else
        assetPath = threeSensorPath;
    end
else
    assetPath = char(options.AssetPath);
end
model = pg_fe_ring_modal_model(string(assetPath), ...
    'ModalDampingRatio',options.ModalDampingRatio);
p = pg_parameters_paper_v3;
p.operating.speedProfile.enabled = false;
assert(p.map.n==18,'The source model must remain exactly 18 DOF.');

if strlength(options.OutputDir)==0
    stamp = datestr(now,'yyyymmdd_HHMMSS');
    outDir = fullfile(projectRoot,'results', ...
        ['fe_causal_path_comparison_' stamp]);
else
    outDir = char(options.OutputDir);
end
if ~exist(outDir,'dir'); mkdir(outDir); end

record = options.RecordDuration;
discard = options.DiscardDuration;
total = record+discard;
fprintf('Causal FE path comparison: fs=%.0f Hz, record=%.6g s, discard=%.6g s\n', ...
    options.fs,record,discard);

%% 2. Solve the two source cases and project PR forces through FE modes
depths = [0 options.CrackDepthMm];
responses = cell(1,numel(depths));
time = [];
first = [];
modalEnergyRecord = [];
faultEventTruth = [];
for ic = 1:numel(depths)
    depth = depths(ic);
    if depth==0
        fault = pg_fault_case('none',0,p);
    else
        fault = pg_root_crack_case(depth,p, ...
            'CrackAngleDeg',options.CrackAngleDeg, ...
            'FaceCoverage',1.0);
    end
    fprintf('  source case q_c=%.3g mm (%d/%d)\n',depth,ic,numel(depths));
    rigid = pg_simulate_lpm_response(total,p,fault,'fs',options.fs);
    motion = pg_motion_state(rigid.time,p);
    % Start the FE modes from the static response to the first moving-load
    % sample.  This removes an artificial zero-state burst while retaining
    % the causal modal response during the record.
    qModal0 = zeros(1,numel(model.frequencyHz));
    for ip0=1:p.model.nPlanet
        theta0 = pg_fe_ring_contact_angle(motion.phi_c(1),p,model,ip0);
        psi0 = pg_fe_ring_input_shape(theta0,model);
        qModal0 = qModal0 + rigid.meshForceRingPlanet(1,ip0).*psi0;
    end
    eta0 = qModal0./(model.omega(:).'.^2);
    causal = pg_apply_fe_ring_transfer_causal(rigid.time, ...
        rigid.meshForceRingPlanet,motion.phi_c(:),p,model, ...
        'InitialModalDisplacement',eta0, ...
        'InitialModalVelocity',zeros(1,numel(eta0)));
    first = round(discard*options.fs)+1;
    last = first+round(record*options.fs)-1;
    if last>numel(rigid.time)
        error('run_causal_fe_ring_path_comparison:Window', ...
            'Requested record exceeds source simulation.');
    end
    if isempty(time)
        time = rigid.time(first:last)-rigid.time(first);
    end
    responses{ic} = causal.acceleration(first:last,:);
    if ic==1
        modalEnergyRecord=causal.modalEnergy(first:last,:);
    else
        faultEventTruth=rigid.faultEventTruth;
        faultEventTruth.eventTimeLocal = faultEventTruth.eventTime - rigid.time(first);
        faultEventTruth.eventStartTimeLocal = faultEventTruth.eventStartTime - rigid.time(first);
        faultEventTruth.eventEndTimeLocal = faultEventTruth.eventEndTime - rigid.time(first);
    end
    clear causal rigid
end

faultIncrement = responses{2}-responses{1};
%% 3. Export signed response and modal-path diagnostics
angles = model.sensorAnglesLocalDeg(:).';
idx0 = find(abs(angles-0)<1e-9,1);
idx90 = find(abs(angles-90)<1e-9,1);
if isempty(idx0)
    error('run_causal_fe_ring_path_comparison:Sensors', ...
        'The comparison figure requires a local 0-degree sensor.');
end
if isempty(idx90)
    idxCompare = find(abs(angles-120)<1e-9,1);
else
    idxCompare = idx90;
end
if isempty(idxCompare)
    error('run_causal_fe_ring_path_comparison:Sensors', ...
        'The comparison figure requires a local 90- or 120-degree sensor.');
end
T = table(time,'VariableNames',{'time_s'});
for ic=1:numel(depths)
    token = sprintf('q%04dum',round(depths(ic)*1000));
    for is=1:numel(angles)
        T.(sprintf('%s_sensor_%03ddeg_acceleration_ms2',token,round(angles(is)))) = ...
            responses{ic}(:,is);
    end
end
for is=1:numel(angles)
    T.(sprintf('q%04dum_minus_healthy_sensor_%03ddeg_acceleration_ms2', ...
        round(options.CrackDepthMm*1000),round(angles(is)))) = faultIncrement(:,is);
end
writetable(T,fullfile(outDir,'causal_fe_path_time_signals.csv'));

% Event-independent path energy report based directly on the FE mode shapes.
E = table(time,'VariableNames',{'time_s'});
for is=1:numel(angles)
    E.(sprintf('sensor_%03ddeg_modal_energy_gain',round(angles(is)))) = ...
        modalEnergyRecord(:,is);
end
writetable(E,fullfile(outDir,'causal_fe_path_modal_energy.csv'));

%% 4. Render comparison figures without changing the response
% A compact event-aligned view: it shows the same source event at all
% sensors, without changing or aligning the waveforms artificially.
fig = figure('Color','w','Visible','off','Position',[100 100 1500 900]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
plot(time,responses{1}(:,idx0)*1e3,'k-','LineWidth',0.7,'DisplayName','healthy 0 deg');
plot(time,responses{2}(:,idx0)*1e3,'b-','LineWidth',0.8,'DisplayName','crack 0 deg');
plot(time,responses{2}(:,idxCompare)*1e3,'r-','LineWidth',0.8, ...
    'DisplayName',sprintf('crack %g deg',angles(idxCompare)));
grid on; box on; xlabel('Time / s'); ylabel('Acceleration / 10^{-3} m s^{-2}');
title(sprintf('Causal FE path response, q_c=%.2f mm',options.CrackDepthMm));
legend('Location','northoutside','NumColumns',3);
nexttile; hold on;
colors=lines(numel(angles));
for is=1:numel(angles)
    plot(time,faultIncrement(:,is)*1e3,'Color',colors(is,:), ...
        'LineWidth',0.8,'DisplayName',sprintf('Delta sensor %g deg',angles(is)));
end
grid on; box on; xlabel('Time / s'); ylabel('Fault increment / 10^{-3} m s^{-2}');
title('Healthy-subtracted causal response (no artificial time alignment)');
legend('Location','northoutside','NumColumns',2);
exportgraphics(fig,fullfile(outDir,'Fig_causal_FE_path_time_response.png'),'Resolution',300);
exportgraphics(fig,fullfile(outDir,'Fig_causal_FE_path_time_response.pdf'),'ContentType','vector');
close(fig);

% Basic event statistics, retaining signs and timing.
stats = table(angles(:),'VariableNames',{'sensor_angle_deg'});
stats.healthy_rms = sqrt(mean(responses{1}.^2,1)).';
stats.fault_rms = sqrt(mean(responses{2}.^2,1)).';
stats.increment_rms = sqrt(mean(faultIncrement.^2,1)).';
stats.increment_peak_abs = max(abs(faultIncrement),[],1).';
writetable(stats,fullfile(outDir,'causal_fe_path_metrics.csv'));

%% 5. Audit fault-event locations and sensor response timing
% Event-aligned timing audit.  No waveform is shifted: each window keeps
% the original local time and the source TVMS event time is only marked.
eventPeakAll = faultEventTruth.eventTimeLocal(:);
eventStartAll = faultEventTruth.eventStartTimeLocal(:);
eventEndAll = faultEventTruth.eventEndTimeLocal(:);
keepEvent = eventStartAll>=0 & eventEndAll<=time(end);
eventTimes = eventPeakAll(keepEvent);
eventStartTimes = eventStartAll(keepEvent);
eventEndTimes = eventEndAll(keepEvent);
eventPlanetIndex = faultEventTruth.planetIndex(keepEvent);
eventCarrierAngleRad = faultEventTruth.carrierAngleRad(keepEvent);
eventContactAngleDeg = zeros(numel(eventTimes),1);
for ie0=1:numel(eventTimes)
    eventContactAngleDeg(ie0)=rad2deg(pg_fe_ring_contact_angle( ...
        eventCarrierAngleRad(ie0),p,model,eventPlanetIndex(ie0)));
end
nEvent = numel(eventTimes);
eventRows = nEvent*numel(angles);
eventNo=zeros(eventRows,1); eventSensor=zeros(eventRows,1);
eventPlanet=zeros(eventRows,1); contactAngle=zeros(eventRows,1);
angularDistance=zeros(eventRows,1);
eventSource=zeros(eventRows,1); eventPeakSource=zeros(eventRows,1);
eventEndSource=zeros(eventRows,1); firstCross=zeros(eventRows,1);
peakTime=zeros(eventRows,1); peakAbs=zeros(eventRows,1); energy5=zeros(eventRows,1);
row=0; halfWindow=0.004;
for ie=1:nEvent
    t0=eventStartTimes(ie); tPeak=eventTimes(ie); tEnd=eventEndTimes(ie);
    ii=find(time>=max(0,t0-halfWindow) & time<=min(time(end),tEnd+halfWindow));
    for is=1:numel(angles)
        row=row+1; x=faultIncrement(ii,is); tt=time(ii);
        [peakAbs(row),ipk]=max(abs(x)); peakTime(row)=tt(ipk)-t0;
        energy=cumsum(x.^2); totalEnergy=energy(end);
        if totalEnergy>0
            ie5=find(energy>=0.05*totalEnergy,1,'first');
            energy5(row)=tt(ie5)-t0;
        else
            energy5(row)=NaN;
        end
        threshold=0.05*peakAbs(row);
        iStart=find(tt>=t0 & abs(x)>=threshold,1,'first');
        if isempty(iStart); firstCross(row)=NaN; else; firstCross(row)=tt(iStart)-t0; end
        eventNo(row)=ie; eventSensor(row)=angles(is); eventSource(row)=t0;
        eventPlanet(row)=eventPlanetIndex(ie); contactAngle(row)=eventContactAngleDeg(ie);
        angularDistance(row)=abs(rad2deg(atan2(sin(deg2rad(contactAngle(row)- ...
            model.sensorAnglesGlobalDeg(is))),cos(deg2rad(contactAngle(row)- ...
            model.sensorAnglesGlobalDeg(is))))));
        eventPeakSource(row)=tPeak; eventEndSource(row)=tEnd;
    end
end
eventMetrics=table(eventNo,eventPlanet,contactAngle,eventSensor,angularDistance, ...
    eventSource,eventPeakSource,eventEndSource, ...
    firstCross,energy5,peakTime,peakAbs, ...
    'VariableNames',{'event_number','active_planet','contact_angle_global_deg', ...
    'sensor_angle_local_deg','contact_sensor_angular_distance_deg', ...
    'source_event_start_time_s', ...
    'source_event_peak_time_s','source_event_end_time_s', ...
    'first_5pct_peak_crossing_lag_s','five_pct_energy_lag_s', ...
    'peak_lag_s','increment_peak_abs_ms2'});
writetable(eventMetrics,fullfile(outDir,'causal_event_arrival_metrics.csv'));

% First complete event, shown without shifting channels, for direct visual
% inspection of the arrival order and modal ringing.
if nEvent>0
    t0=eventStartTimes(1); tPeak=eventTimes(1); tEnd=eventEndTimes(1);
    ii=find(time>=max(0,t0-halfWindow) & time<=min(time(end),tEnd+halfWindow));
    fig2=figure('Color','w','Visible','off','Position',[120 120 1400 780]);
    hold on; colors=lines(numel(angles));
    for is=1:numel(angles)
        plot(time(ii)-t0,faultIncrement(ii,is)*1e3,'Color',colors(is,:), ...
            'LineWidth',0.9,'DisplayName',sprintf('%g deg',angles(is)));
    end
    xline(0,'k--','LineWidth',0.8,'DisplayName','fault contact starts');
    xline(tPeak-t0,'k:','LineWidth',0.8,'DisplayName','maximum stiffness loss');
    grid on; box on; xlabel('Time from source event / s');
    ylabel('Fault increment / 10^{-3} m s^{-2}');
    title(sprintf('Unshifted first fault event, q_c=%.2f mm',options.CrackDepthMm));
    legend('Location','northeast');
    exportgraphics(fig2,fullfile(outDir,'Fig_causal_FE_path_first_event_unshifted.png'),'Resolution',300);
    exportgraphics(fig2,fullfile(outDir,'Fig_causal_FE_path_first_event_unshifted.pdf'),'ContentType','vector');
    close(fig2);
end

%% 6. Write provenance and return the computed signals
fid=fopen(fullfile(outDir,'README.txt'),'w');
fprintf(fid,'Causal moving-contact FE modal path comparison\n');
fprintf(fid,'Source: existing 18-DOF LPM/Newmark model; no source parameters changed.\n');
fprintf(fid,'FE asset: %s\n',assetPath);
fprintf(fid,'Modes: %d, %.9g..%.9g Hz; damping ratio=%.9g provisional.\n', ...
    numel(model.frequencyHz),min(model.frequencyHz),max(model.frequencyHz), ...
    options.ModalDampingRatio);
fprintf(fid,'Integrator: exact zero-order-hold modal state transition; initial state=static.\n');
fprintf(fid,'fs=%.9g Hz, discard=%.9g s, record=%.9g s.\n',options.fs,discard,record);
fprintf(fid,'The moving contact shape and sensor mode shapes provide all spatial attenuation; no distance gain, manual delay, smoothing, or channel scaling was applied.\n');
fprintf(fid,'The current FE asset is a constrained ring-box modal asset; SI unit-mass normalization is unverified, so absolute acceleration and FRF magnitudes are exploratory.\n');
fprintf(fid,'For full 15--45 kHz causal timing, fs=524288 Hz is a diagnostic setting; verify convergence at higher fs before quantitative claims.\n');
fprintf(fid,'Event timing columns are calculated relative to the source TVMS event without channel alignment; maximum ringing peaks are not first-arrival times.\n');
fclose(fid);

result=struct('outputDir',outDir,'time',time,'depths',depths, ...
    'angles',angles,'healthy',responses{1},'fault',responses{2}, ...
    'faultIncrement',faultIncrement,'metrics',stats);
fprintf('Causal FE path results written to %s\n',outDir);
end
