function result = run_tvms
%RUN_TVMS Healthy TVMS entry point.  It refuses unresolved geometry.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot,'src'),fullfile(projectRoot,'configs'), ...
    fullfile(projectRoot,'tvms'));
p = pg_parameters;
c = tvms_config(p);
audit = tvms_geometry_audit(c);
outDir = fullfile(projectRoot,'results','tvms');
if ~isfolder(outDir), mkdir(outDir); end
fid = fopen(fullfile(outDir,'tvms_summary.txt'),'w');
fprintf(fid,'Healthy potential-energy TVMS preflight\n\n');
fprintf(fid,'Zs=%d; Zp=%d; Zr=%d; m=%.9g m; alpha=%.9g rad\n',c.zs,c.zp,c.zr,c.module,c.alphaNominal);
fprintf(fid,'zero-modification reference Zr=Zs+2Zp=%d\n',audit.standardRingTeeth);
fprintf(fid,'centre modifications: ySP=%+.9g; yPR=%+.9g; centre-modified Zr=%.9g; residual=%+.9g\n', ...
    audit.centerModificationSP,audit.centerModificationPR, ...
    audit.centerModifiedRingTeeth,audit.toothClosureResidual);
fprintf(fid,'SP profile shift: supplied sum=%+.9g; centre-implied sum=%+.9g; residual=%+.9g\n', ...
    audit.suppliedProfileShiftSumSP,audit.impliedProfileShiftSumSP, ...
    audit.profileShiftResidualSP);
fprintf(fid,'PR profile shift: supplied difference=%+.9g; centre-implied difference=%+.9g; residual=%+.9g\n', ...
    audit.suppliedProfileShiftDifferencePR,audit.impliedProfileShiftDifferencePR, ...
    audit.profileShiftResidualPR);
fprintf(fid,'external center residual=%.9g m; internal center residual=%.9g m\n',audit.externalCenterResidual,audit.internalCenterResidual);
fprintf(fid,'%s\n',audit.warning);
if ~isempty(audit.missing), fprintf(fid,'TODO required before computation: %s\n',strjoin(audit.missing,', ')); end
fclose(fid);
if ~audit.ready
    error('run_tvms:UnresolvedGeometry',['TVMS not computed. %s Required values are recorded in configs/tvms_config.m and results/tvms/tvms_summary.txt.'],audit.warning);
end
error('run_tvms:NotImplemented','Potential-energy evaluator must not run until the reviewed geometric TODO values are supplied.');
end
