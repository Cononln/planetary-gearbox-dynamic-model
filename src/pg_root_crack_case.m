function fault = pg_root_crack_case(depthMm,p,options)
%PG_ROOT_CRACK_CASE Energy-based single sun-tooth root-crack definition.
% depthMm is the crack extension q in millimetres. The default 65-deg
% crack angle and 0.5/1.5/2.5-mm study levels follow Shen (2023).
%
% SpeedRippleFraction is the amplitude a_m/n_0 in
%   n(t)=n_0+a_m*cos(2*pi*f_s*t+phase).
% Shen did not report a_m, so the default is zero and the optional imposed
% phase modulation remains auditable rather than being silently fitted.

arguments
    depthMm (1,1) double {mustBeNonnegative}
    p struct = pg_parameters
    options.CrackAngleDeg (1,1) double {mustBeGreaterThan(options.CrackAngleDeg,0),mustBeLessThan(options.CrackAngleDeg,90)} = p.tvms.rootCrack.angleDeg
    options.FaceCoverage (1,1) double {mustBeGreaterThanOrEqual(options.FaceCoverage,0),mustBeLessThanOrEqual(options.FaceCoverage,1)} = p.tvms.rootCrack.faceCoverage
    options.SpeedRippleFraction (1,1) double {mustBeNonnegative} = 0
    options.SpeedRipplePhase (1,1) double = 0
end

fault.type = 'sun';
fault.damageModel = 'root_crack_energy';
fault.severity = depthMm;
fault.tooth = p.tvms.faultTooth;
fault.targetPlanet = p.tvms.targetPlanet;
fault.crackDepth = depthMm*1e-3;
fault.crackDepthMm = depthMm;
fault.depthToModuleRatio = fault.crackDepth/p.gear.module;
fault.crackAngleDeg = options.CrackAngleDeg;
fault.crackAngle = deg2rad(options.CrackAngleDeg);
fault.faceCoverage = options.FaceCoverage;
fault.fullFaceWidth = p.gear.faceWidthSunPlanet;
fault.crackedFaceWidth = options.FaceCoverage*fault.fullFaceWidth;
fault.speedRippleFraction = options.SpeedRippleFraction;
fault.speedRipplePhase = options.SpeedRipplePhase;
fault.depthDefinition = 'tooth-root crack extension q';

[contactFraction,ratio,energy] = pg_root_crack_energy_lookup(p,fault);
fault.energyLookup.contactFraction = contactFraction;
fault.energyLookup.stiffnessRatio = ratio;
fault.energyLookup.detail = energy;
fault.nearSeveredNumericalLimit = min(ratio) <= ...
    p.tvms.rootCrack.minLigamentFraction+1e-12;
end
