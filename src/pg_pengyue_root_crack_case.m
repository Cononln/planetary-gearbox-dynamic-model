function fault = pg_pengyue_root_crack_case(depthMm,angleDeg,p)
%PG_PENGYUE_ROOT_CRACK_CASE Straight sun-root crack used in Peng Yue Fig.3-14.

arguments
    depthMm (1,1) double {mustBePositive}
    angleDeg (1,1) double {mustBeGreaterThan(angleDeg,0),mustBeLessThan(angleDeg,90)}
    p struct = pg_parameters_pengyue_fig314
end

fault.type = 'sun';
fault.damageModel = 'root_crack_energy';
fault.severity = depthMm;
fault.tooth = p.tvms.faultTooth;
fault.targetPlanet = p.tvms.targetPlanet;
fault.crackDepth = depthMm*1e-3;
fault.crackDepthMm = depthMm;
fault.depthToModuleRatio = fault.crackDepth/p.gear.module;
fault.crackAngleDeg = angleDeg;
fault.crackAngle = deg2rad(angleDeg);
fault.faceCoverage = 1.0;
fault.fullFaceWidth = p.gear.faceWidthSunPlanet;
fault.crackedFaceWidth = fault.fullFaceWidth;
fault.speedRippleFraction = 0;
fault.speedRipplePhase = 0;
fault.depthDefinition = 'straight tooth-root crack extension q';

[contactFraction,ratio,energy] = pg_pengyue_root_crack_energy_lookup(p,fault);
fault.energyLookup.contactFraction = contactFraction;
fault.energyLookup.stiffnessRatio = ratio;
fault.energyLookup.detail = energy;
fault.nearSeveredNumericalLimit = min(ratio) <= ...
    p.tvms.rootCrack.minLigamentFraction+1e-12;
end
