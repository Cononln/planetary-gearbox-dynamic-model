function report = diagnose_pengyue_foundation_q3(options)
%DIAGNOSE_PENGYUE_FOUNDATION_Q3 TVMS-only foundation sensitivity audit.
% Contact-ratio and active-pair logic are unchanged. The 7e-9 m/N case is
% a temporary constant-compliance sensitivity, not a final Chaari geometry.
arguments
    options.OutputDir (1,1) string = ""
    options.Samples (1,1) double {mustBeInteger,mustBeGreaterThan(options.Samples,1000)} = 4001
    options.FoundationDiagnosticCompliance (1,1) double {mustBeNonnegative} = 7e-9
end
setup_paths;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
if strlength(options.OutputDir)==0
    outDir = fullfile(projectRoot,'results','pengyue_tvms_foundation_diagnostic_q3');
else
    outDir = char(options.OutputDir);
end
if ~isfolder(outDir); mkdir(outDir); end
p = pg_parameters_pengyue_fig314;
assert(p.map.n==18,'Expected 18 DOF.');
fH = struct('type','none','tooth',1,'crackDepth',0,'crackAngle',0);
fF = struct('type','sun','tooth',1,'crackDepth',3e-3,'crackAngle',deg2rad(45),'crackAngleDeg',45);
t = linspace(0,1/p.kin.f_mesh,options.Samples).';
cases = {'healthy',fH;'q3_gamma45',fF};
foundation = [0,options.FoundationDiagnosticCompliance];
rows = table(); data = struct();
for jf=1:numel(foundation)
    c = foundation(jf);
    if c==0, fModel="none"; tag='foundation_none'; else, fModel="constant"; tag='foundation_constant'; end
    for jc=1:size(cases,1)
        name=cases{jc,1}; fault=cases{jc,2};
        mdl=pg_build_pengyue_sp_tvms_model(p,fault,'FoundationModel',fModel, ...
            'FoundationComplianceSun',c,'FoundationCompliancePlanet',c, ...
            'LookupPoints',1601,'IntegrationPoints',1200);
        [k,d]=pg_eval_pengyue_sp_tvms(t,p,fault,mdl); y=k(:,1); np=d.activePairCount(:,1);
        one=np==1; two=np==2;
        row=table(string(tag),string(name),c,min(y),max(y),mean(y),mean(y(one)),mean(y(two)),mean(one),mean(two),mdl.contactRatio, ...
            'VariableNames',{'foundation_model','case','foundation_compliance_m_per_N','k_min_N_per_m','k_max_N_per_m','k_mean_N_per_m','single_mean_N_per_m','double_mean_N_per_m','single_fraction','double_fraction','contact_ratio'});
        rows=[rows;row]; %#ok<AGROW>
        data.(tag).(name).phase_norm=t*p.kin.f_mesh; data.(tag).(name).k_sp1_N_per_m=y; data.(tag).(name).active_pair_count=np; data.(tag).(name).model=mdl;
        fprintf('%s | %s | min %.6e max %.6e mean %.6e | single %.6e double %.6e\n',tag,name,min(y),max(y),mean(y),mean(y(one)),mean(y(two)));
    end
end
csvPath=fullfile(outDir,'foundation_q3_tvms_stats.csv'); matPath=fullfile(outDir,'foundation_q3_tvms.mat'); txtPath=fullfile(outDir,'foundation_q3_tvms_diagnostic.txt');
writetable(rows,csvPath); save(matPath,'rows','data','p','-v7.3');
fid=fopen(txtPath,'w'); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Peng Yue SP TVMS foundation diagnostic\nzs=%d,zp=%d,zr=%d,m=%.9g,alpha=%.9g rad\n',p.gear.zs,p.gear.zp,p.gear.zr,p.gear.module,p.gear.alpha_nominal);
fprintf(fid,'ns=%.9g rpm,fm=%.12g Hz,q=3 mm,gamma=45 deg\ncontact ratio=%.12g; active-pair logic unchanged\n',p.operating.sunRpm,p.kin.f_mesh,rows.contact_ratio(1));
fprintf(fid,'The 7e-9 m/N rows are temporary constant-compliance sensitivity only.\n');
fprintf(fid,'Final Chaari run requires explicit r_int and theta_f for both gears; current source status is MISSING.\n\n');
for i=1:height(rows)
    fprintf(fid,'%s,%s,c=%.9e,min=%.9e,max=%.9e,mean=%.9e,single=%.9e,double=%.9e,singleFrac=%.9f,doubleFrac=%.9f\n',rows.foundation_model(i),rows.case(i),rows.foundation_compliance_m_per_N(i),rows.k_min_N_per_m(i),rows.k_max_N_per_m(i),rows.k_mean_N_per_m(i),rows.single_mean_N_per_m(i),rows.double_mean_N_per_m(i),rows.single_fraction(i),rows.double_fraction(i));
end
fprintf(fid,'\nCSV: %s\nMAT: %s\nTXT: %s\n',csvPath,matPath,txtPath);
report=struct('rows',rows,'csv',csvPath,'mat',matPath,'txt',txtPath);
fprintf('Wrote %s\nWrote %s\nWrote %s\n',csvPath,matPath,txtPath);
end
