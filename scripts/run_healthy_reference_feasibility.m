function result=run_healthy_reference_feasibility(commonMode)
%RUN_HEALTHY_REFERENCE_FEASIBILITY Fixed-protocol healthy-to-fault pilot.
% Numerical analysis only; figure generation is a separate Python script.
setup_paths;
if nargin<1, commonMode="incremental"; end
commonMode=string(commonMode);
assert(any(commonMode==["incremental","calibration_anchored","nominal"]));
root=fileparts(fileparts(mfilename('fullpath')));
outDir=fullfile(root,'results','healthy_reference_feasibility_20260905');
if commonMode~="incremental", outDir=outDir+"_exploratory_"+commonMode; end
if isfile(fullfile(outDir,'pilot_result.mat'))
    error('Output exists. Preserve it; select a new run directory for a new experiment.');
end
if ~isfolder(outDir), mkdir(outDir); end
train=[3 15]; validation=[17 25]; blocks=[27 37;38 48;49 59];
files={fullfile(root,'data_local','bl','healthy_record.MAT'), ...
    fullfile(root,'data_local','planet','pf50_record.MAT')};
conditions=["BL","PF50"];
fs=5120; zr=84; h=1:4;
records=cell(1,2);
for c=1:2
    fprintf('Loading / finite-support preprocessing %s...\n',conditions(c));
    records{c}=preprocess(files{c},fs,h,zr,train,commonMode);
end
healthy=records{1}; fault=records{2};
it=interval(healthy.t,train); iv=interval(healthy.t,validation);
zz=healthy.z(it,:,:); tt=healthy.common.carrierPhase(it);
scale=max(median(abs(zz),1),eps);
harmonicScale=sqrt(mean(mean(abs(zz).^2,1),2));
spec=struct('name',{},'orders',{},'penalty',{});
spec(1)=struct('name',"constant",'orders',0,'penalty',1);
spec(2)=struct('name',"strict3k",'orders',-12:3:12,'penalty',ones(1,9));
for k=[3 6 12]
    spec(end+1)=struct('name',"free"+k,'orders',-k:k,'penalty',ones(1,2*k+1)); %#ok<AGROW>
end
for p=[10 100]
    orders=-12:12; penalty=ones(size(orders)); penalty(mod(orders,3)~=0)=p;
    spec(end+1)=struct('name',"soft3k_"+p,'orders',orders,'penalty',penalty); %#ok<AGROW>
end
subsets={1:3,[1 2],[1 3],[2 3]}; models={}; candidateRows={};
fprintf('Selecting models exclusively on healthy validation...\n');
for s=1:numel(subsets)
    for k=1:numel(spec)
        model=pg_frozen_path_fit(zz,tt,subsets{s},spec(k).orders,spec(k).penalty,.002);
        idx=find(iv); idx=idx(1:8:end); y=pg_frozen_path_predict(model,healthy.z(idx,:,:),healthy.common.carrierPhase(idx),0);
        base=pg_phase_resultant(healthy.z(idx,:,:),model.subset,scale);
        corrected=pg_phase_resultant(y,model.subset,scale);
        gain=zeros(1,2);
        for half=1:2
            use=(healthy.t(idx)>=validation(1)+(half-1)*diff(validation)/2) & ...
                (healthy.t(idx)<validation(1)+half*diff(validation)/2);
            gain(half)=pg_phase_resultant(y(use,:,:),model.subset,scale)- ...
                pg_phase_resultant(healthy.z(idx(use),:,:),model.subset,scale);
        end
        model.specName=spec(k).name;
        models{end+1}=model; %#ok<AGROW>
        candidateRows(end+1,:)={numel(models),s,spec(k).name,join(string(model.subset),"+"), ...
            base,corrected,corrected-base,min(gain), ...
            corrected-base>=.015 && min(gain)>=-.01}; %#ok<AGROW>
    end
end
candidates=cell2table(candidateRows,'VariableNames',{'id','subset_id','specification','sensors', ...
    'base_coherence','corrected_coherence','gain','worst_half_gain','admissible'});
writetable(candidates,fullfile(outDir,'healthy_validation_candidates.csv'));
identity=pg_frozen_path_fit(zz,tt,1:3,[],[],.002); identity.specName="identity";
selected=selectModel(candidates,models,identity,false);
allModel=selectModel(candidates,models,identity,true);
constant=pg_frozen_path_fit(zz,tt,selected.subset,0,1,.002); constant.specName="constant";
fprintf('Healthy selection: %s sensors %s. All-channel: %s\n',selected.specName,mat2str(selected.subset),allModel.specName);
faultCal=interval(fault.t,train);
offsetSelected=pg_frozen_path_offset(selected,fault.z(faultCal,:,:),fault.common.carrierPhase(faultCal));
offsetAll=pg_frozen_path_offset(allModel,fault.z(faultCal,:,:),fault.common.carrierPhase(faultCal));
localFault=pg_frozen_path_fit(fault.z(faultCal,:,:),fault.common.carrierPhase(faultCal), ...
    selected.subset,selected.orders,selected.penalty,.002);
localFault.specName=selected.specName;
frozenBefore=selected;
methodNames=["common_all","common_selected","healthy_constant","healthy_frozen_all", ...
    "healthy_frozen_selected","within_record","sensor1","sensor2","sensor3"];
metricRows={}; spectrumRows={}; angleRows={}; encoderRows={};
for c=1:2
    r=records{c};
    % Encoder use begins here, after all vibration-only inference / selection.
    valid=r.t>=3 & r.t<59;
    encoderOffset=median(r.common.carrierPhase(valid)-r.encoderA(valid));
    nominal=2*pi*r.fmesh/zr*r.t;
    nominalOffset=median(nominal(valid)-r.encoderA(valid));
    encAB=r.encoderA-r.encoderB; encAB=encAB-median(encAB(valid));
    speedEncoder=movmean(gradient(r.encoderA)*fs*60/(2*pi),round(.1*fs));
    speedEstimate=r.common.frequencyHz/zr*60;
    for b=1:size(blocks,1)
        use=interval(r.t,blocks(b,:)); z=r.z(use,:,:); theta=r.common.carrierPhase(use);
        e=(theta-r.encoderA(use)-encoderOffset)*180/pi;
        en=(nominal(use)-r.encoderA(use)-nominalOffset)*180/pi;
        encoderRows(end+1,:)={conditions(c),b,r.fmesh,rms(e),rms(en), ...
            rms(speedEstimate(use)-speedEncoder(use)),rms(encAB(use))*180/pi, ...
            mean(speedEncoder(use)),mean(speedEstimate(use)),r.encoderMeta.edgeRateA,r.encoderMeta.edgeRateB}; %#ok<AGROW>
        ix=find(use); ix=ix(1:52:end);
        angleRows{end+1}=table(repmat(conditions(c),numel(ix),1),repmat(b,numel(ix),1),r.t(ix), ...
            (r.common.carrierPhase(ix)-r.encoderA(ix)-encoderOffset)*180/pi, ...
            (nominal(ix)-r.encoderA(ix)-nominalOffset)*180/pi,speedEstimate(ix),speedEncoder(ix), ...
            'VariableNames',{'condition','block','time_s','estimated_error_deg','nominal_error_deg','estimated_rpm','encoder_rpm'}); %#ok<AGROW>
        rawDiff=envelopes(z,1:3,harmonicScale,"differential");
        [diffScore,diffSpectrum]=scoreEnvelope(rawDiff,theta,zr);
        spectrumRows{end+1}=spectrumTable(diffSpectrum,conditions(c),b,"retained_raw_differential"); %#ok<AGROW>
        for m=1:numel(methodNames)
            subset=1:3; y=z;
            switch m
                case 2, subset=selected.subset;
                case 3
                    subset=selected.subset; y=pg_frozen_path_predict(constant,z,theta,0);
                case 4
                    offset=0; if c==2, offset=offsetAll; end
                    y=pg_frozen_path_predict(allModel,z,theta,offset);
                case 5
                    subset=selected.subset; offset=0; if c==2, offset=offsetSelected; end
                    y=pg_frozen_path_predict(selected,z,theta,offset);
                case 6
                    subset=selected.subset; model=selected; if c==2, model=localFault; end
                    y=pg_frozen_path_predict(model,z,theta,0);
                case {7,8,9}, subset=m-6;
            end
            [coh,unweighted]=pg_phase_resultant(y,subset,scale);
            envelope=envelopes(y,subset,harmonicScale,"common");
            [score,spectrum]=scoreEnvelope(envelope,theta,zr);
            ampError=max(abs(abs(y(:))-abs(z(:))))/max(abs(z(:)));
            metricRows(end+1,:)={conditions(c),b,methodNames(m),join(string(subset),"+"), ...
                coh,unweighted,score(1),score(2),diffScore(1),diffScore(2),ampError}; %#ok<AGROW>
            spectrumRows{end+1}=spectrumTable(spectrum,conditions(c),b,methodNames(m)); %#ok<AGROW>
        end
    end
end
assert(isequaln(selected,frozenBefore),'Frozen coefficients changed.');
metrics=cell2table(metricRows,'VariableNames',{'condition','block','method','sensors','coherence', ...
    'unweighted_coherence','fp_family_db','twofp_family_db','retained_diff_fp_db','retained_diff_twofp_db','amplitude_error'});
encoderMetrics=cell2table(encoderRows,'VariableNames',{'condition','block','mesh_hz','carrier_angle_rmse_deg', ...
    'nominal_angle_rmse_deg','carrier_speed_rmse_rpm','encoder_AB_angle_rmse_deg', ...
    'encoder_mean_rpm','estimated_mean_rpm','encoder_A_edges_per_s','encoder_B_edges_per_s'});
writetable(metrics,fullfile(outDir,'heldout_metrics.csv'));
writetable(encoderMetrics,fullfile(outDir,'encoder_validation.csv'));
writetable(vertcat(spectrumRows{:}),fullfile(outDir,'envelope_spectra.csv'));
writetable(vertcat(angleRows{:}),fullfile(outDir,'angle_trace.csv'));
result.metrics=metrics; result.encoderMetrics=encoderMetrics; result.candidates=candidates;
result.selected=selected; result.allModel=allModel; result.constant=constant; result.localFault=localFault;
result.offsetSelectedRad=offsetSelected; result.offsetAllRad=offsetAll;
result.blocks=blocks; result.train=train; result.validation=validation;
result.fs=fs; result.zr=zr; result.harmonics=h; result.harmonicScale=harmonicScale;
result.file=files; result.encoderMeta={records{1}.encoderMeta,records{2}.encoderMeta};
result.sourceRate=[records{1}.sourceFs,records{2}.sourceFs];
result.commonMode=commonMode;
result.exploratory=commonMode~="incremental";
result.calibrationSlopeBiasHz=[records{1}.calibrationSlopeBiasHz,records{2}.calibrationSlopeBiasHz];
result.modelFrozen=isequaln(selected,frozenBefore);
result.limitations='One record per condition; temporal blocks are not independent experimental replicates. Effective relative phase is not a uniquely identified structural path.';
save(fullfile(outDir,'pilot_result.mat'),'result','-v7.3');
fprintf('DONE. Selected offset %.3f deg, all-channel offset %.3f deg.\n',rad2deg(offsetSelected),rad2deg(offsetAll));
disp(metrics); disp(encoderMetrics);
end

function model=selectModel(candidates,models,identity,onlyAll)
ok=candidates.admissible;
if onlyAll, ok=ok & candidates.subset_id==1; end
if ~any(ok), model=identity; return; end
idx=find(ok); [~,best]=max(candidates.corrected_coherence(idx)); model=models{idx(best)};
end

function out=preprocess(file,fs,h,zr,calibration,commonMode)
a=load(file); assert(isfield(a,'Data') && size(a.Data,2)==7);
sourceFs=double(a.SampleRate); assert(size(a.Data,1)/sourceFs>=59.5);
x=resample(double(a.Data(:,1:3)),fs,sourceFs);
t=(0:size(x,1)-1)'/fs;
x=x-mean(x(interval(t,calibration),:),1);
% Pulse decoding is stored for later validation, never used for inference.
[encA,metaA]=decodePulse(a.Data(:,5),sourceFs,t);
[encB,metaB]=decodePulse(a.Data(:,6),sourceFs,t);
clear a
use=interval(t,calibration); xx=x(use,:); n=size(xx,1); w=hann(n);
amp=sqrt(mean(abs(fft(xx.*w)).^2,2)); f=(0:n-1)'*fs/n;
search=find(f>=145 & f<=190); [~,k]=max(amp(search)); i=search(k);
v=log(max(amp(i-1:i+1),eps)); d=.5*(v(1)-v(3))/(v(1)-2*v(2)+v(3));
fmesh=f(i)+max(-1,min(1,d))*fs/n;
N=2*round(.25*fs); delay=N/2; z=complex(zeros(size(x,1),3,numel(h)));
for k=1:numel(h)
    bw=30+4*(h(k)-1); fir=fir1(N,bw/(fs/2),kaiser(N+1,6));
    osc=exp(1i*2*pi*h(k)*fmesh*t);
    base=x.*conj(osc); y=fftfilt(fir,[base;zeros(delay,3)]);
    z(:,:,k)=2*y(delay+1:delay+size(x,1),:).*osc;
end
common=pg_common_phase_frozen(z,fs,h,fmesh,zr,median(abs(z(use,:,:)),1));
nominal=2*pi*fmesh*t;
slope=polyfit(t(use),common.meshPhase(use)-nominal(use),1);
if commonMode=="calibration_anchored"
    common.meshPhase=common.meshPhase-slope(1)*t;
    common.frequencyHz=common.frequencyHz-slope(1)/(2*pi);
elseif commonMode=="nominal"
    common.meshPhase=nominal;
    common.frequencyHz=ones(size(t))*fmesh;
end
common.carrierPhase=common.meshPhase/zr;
common.monotonic=all(diff(common.carrierPhase)>0);
assert(common.monotonic,'Carrier phase lost monotonicity.');
out=struct('z',z,'t',t,'fs',fs,'sourceFs',sourceFs,'fmesh',fmesh,'common',common, ...
    'calibrationSlopeBiasHz',slope(1)/(2*pi), ...
    'encoderA',encA,'encoderB',encB,'encoderMeta',struct('A',metaA,'B',metaB, ...
    'edgeRateA',metaA.edgeRate,'edgeRateB',metaB.edgeRate));
fprintf('  fs=%g -> %g Hz, mesh=%.6f Hz, encoder A/B=%.3f/%.3f edges/s\n',sourceFs,fs,fmesh,metaA.edgeRate,metaB.edgeRate);
end

function [angle,meta]=decodePulse(y,fs,t)
% Two thresholds suppress chatter; use the high-threshold rising crossing.
lo=prctile(y(1:31:end),1); hi=prctile(y(1:31:end),99);
up=lo+.65*(hi-lo); down=lo+.35*(hi-lo);
rise=find(y(1:end-1)<up & y(2:end)>=up)+1;
fall=find(y(1:end-1)>down & y(2:end)<=down)+1;
events=[rise,ones(size(rise));fall,zeros(size(fall))]; events=sortrows(events,1);
keep=events(:,2)==1 & [true;events(1:end-1,2)==0]; edge=events(keep,1);
edge=edge([true;diff(edge)>.00015*fs]);
te=(edge-1)/fs; ae=(0:numel(edge)-1)'*2*pi/1024;
angle=interp1(te,ae,t,'linear','extrap');
meta=struct('low',lo,'high',hi,'edgeCount',numel(edge),'edgeRate',(numel(edge)-1)/(te(end)-te(1)), ...
    'medianIntervalSec',median(diff(te)),'badIntervalFraction',mean(abs(diff(te)/median(diff(te))-1)>.35));
end

function use=interval(t,lim)
use=t>=lim(1) & t<lim(2);
end

function envelope=envelopes(z,subset,scale,branch)
zz=z(:,subset,:);
if branch=="common", a=abs(mean(zz,2));
else, a=sqrt(mean(abs(zz-mean(zz,2)).^2,2)); end
envelope=sum(a./max(scale,eps),3); envelope=envelope(:);
end

function [scores,out]=scoreEnvelope(envelope,theta,zr)
turn=theta/(2*pi); tq=(ceil(turn(1)*256):floor(turn(end)*256))'/256;
e=interp1(turn,envelope,tq,'linear'); e=e-mean(e); n=numel(e); w=hann(n);
amp=2*abs(fft(e.*w))/sum(w); order=(0:n-1)'*256/n;
fp=zr/31; lines=unique([fp*(1:4),2*fp*(1:4)]);
bg=order>=1 & order<=30;
for f=lines, bg=bg & abs(order-f)>.45; end
background=max(median(amp(bg)),eps); scores=zeros(1,2);
for family=1:2
    peaks=zeros(1,4);
    for k=1:4, peaks(k)=max(amp(abs(order-family*fp*k)<=.18)); end
    scores(family)=20*log10(max(rms(peaks),eps)/background);
end
keep=order<=30;
out=table(order(keep),amp(keep),20*log10(max(amp(keep),eps)/background), ...
    'VariableNames',{'carrier_order','amplitude','background_relative_db'});
end

function tbl=spectrumTable(s,condition,block,method)
n=height(s); tbl=[table(repmat(condition,n,1),repmat(block,n,1),repmat(method,n,1), ...
    'VariableNames',{'condition','block','method'}),s];
end
