function p = pg_parameters_pengyue_fig314
%PG_PARAMETERS_PENGYUE_FIG314 Parameters for Peng Yue (2025) Fig. 3-14.
% This is an isolated reproduction configuration. It does not modify the
% repository's current 21/31/84 experimental model.
%
% Values marked "reported" are taken directly from Table 2-1 and Sec. 3.4.1
% of Peng Yue, Research on Fault Dynamics Model and Characteristics of
% Planetary Gears, NCEPU, 2025. Values under p.reproduction.assumption are
% necessary numerical/model details not reported in the thesis and must not
% be presented as source-paper facts.

p = pg_parameters;

p.model.name = 'Peng Yue 2025 Fig. 3-14 reproduction';
p.model.nPlanet = 3;
p.model.fs = 10240; % assumption: same rate as the later bench acquisition
p.model.phi_c0 = 0;

% Reported geometry / Table 2-1.
p.gear.zs = 17;
p.gear.zp = 37;
p.gear.zr = 91;
p.gear.module = 5e-3;
p.gear.alpha_nominal = deg2rad(20);
p.gear.planetPhase0 = 2*pi*(0:p.model.nPlanet-1)/p.model.nPlanet;

p.gear.pitchRadiusSun = p.gear.module*p.gear.zs/2;
p.gear.pitchRadiusPlanet = p.gear.module*p.gear.zp/2;
p.gear.pitchRadiusRing = p.gear.module*p.gear.zr/2;
% Table 2-1 reports the following base-circle diameters explicitly.
p.gear.baseRadiusSun = 79.8793e-3/2;
p.gear.baseRadiusPlanet = 173.8431e-3/2;
p.gear.baseRadiusRing = 427.5601e-3/2;
% Physical planet-centre radius from the standard geometry. The table also
% reports 253.7 mm for the carrier equivalent base diameter; keep it
% separately for converting the reported torsional support parameter.
p.gear.carrierRadius = p.gear.pitchRadiusSun+p.gear.pitchRadiusPlanet;
p.gear.carrierTorsionalRadius = 253.7e-3/2;

% Standard full-depth involute radii. The thesis gives m and z, but does not
% tabulate addendum/root diameters. These values are therefore reproducible
% geometry assumptions, not directly reported numbers.
p.gear.addendumRadiusSun = p.gear.module*(p.gear.zs+2)/2;
p.gear.rootRadiusSun = p.gear.module*(p.gear.zs-2.5)/2;
p.gear.addendumRadiusPlanet = p.gear.module*(p.gear.zp+2)/2;
p.gear.rootRadiusPlanet = p.gear.module*(p.gear.zp-2.5)/2;

% Chaari/Sainsot fillet-foundation geometry is not reported in the source
% thesis (the inner/body radius and root half-angle are needed in addition
% to module and tooth count).  Keep the fields explicit and missing rather
% than inferring a solid-body radius or silently fitting a compliance.
p.gear.foundation.sun.innerRadius = NaN;
p.gear.foundation.sun.rootHalfAngle = NaN;
p.gear.foundation.planet.innerRadius = NaN;
p.gear.foundation.planet.rootHalfAngle = NaN;
p.gear.foundation.ring.rimOuterRadius = NaN;
p.gear.foundation.ring.rootHalfAngle = NaN;
p.gear.foundation.ring.filletGeometry = 'MISSING';
p.gear.foundation.ufDefinition = 'explicit Chaari Fig. 2 geometry required';
p.gear.addendumRadiusRing = p.gear.module*(p.gear.zr-2)/2;
p.gear.rootRadiusRing = p.gear.module*(p.gear.zr+2.5)/2;

% The thesis does not report tooth face width or E/nu in Table 2-1. They are
% explicit assumptions here so that they can be changed without touching
% the reproduction logic.
p.gear.faceWidthSunPlanet = 20e-3;
p.gear.faceWidthRing = 20e-3;
p.material.gear.E = 206e9;
p.material.gear.nu = 0.30;
p.material.gear.G = p.material.gear.E/(2*(1+p.material.gear.nu));

% Table 2-1 masses and inertias. Convert kg*mm^2 to kg*m^2.
p.body.sun.mass = 20.29;
p.body.planet.mass = 3.91;
p.body.ring.mass = 59.88;
p.body.carrier.mass = 57.314;
p.body.sun.J = 27049.57e-6;
p.body.planet.J = 22368.33e-6;
p.body.ring.J = 3853845.79e-6;
p.body.carrier.J = 776910.239e-6;

% Reported operating point.
p.operating.sunRpm = 800;
p.operating.encoderPpr = 1024;
p.kin = pg_kinematics(p.gear.zs,p.gear.zp,p.gear.zr, ...
    p.operating.sunRpm,p.operating.encoderPpr);
p.operating.outputTorque = 50; % N m, reported Tout
% Tin is not stated in Sec. 3.4.1. Use ideal steady power balance only as an
% initial-load assumption; this is logged below and can be overridden.
p.operating.inputTorque = p.operating.outputTorque* ...
    abs(p.kin.omega_c/p.kin.omega_s);
p.operating.torqueRippleFraction = 0;
p.operating.torqueRippleFrequencyHz = p.kin.f_s;

% Table 2-1 support stiffness values are N/um. Convert to N/m.
kSun = 10e6;
kPlanet = 58e6;
kRing = 1000e6;
kCarrier = 195e6;
zetaSupport = 0.03; % reported dimensionless support damping coefficient
p.support.sun.kxy = kSun;
p.support.ring.kxy = kRing;
p.support.carrier.kxy = kCarrier;
p.support.planetPin.kxy = kPlanet;
p.support.sun.cxy = supportDamping(kSun,p.body.sun.mass,zetaSupport);
p.support.ring.cxy = supportDamping(kRing,p.body.ring.mass,zetaSupport);
p.support.carrier.cxy = supportDamping(kCarrier,p.body.carrier.mass,zetaSupport);
p.support.planetPin.cxy = supportDamping(kPlanet,p.body.planet.mass,zetaSupport);

% Table 2-1 torsional support stiffness is also listed in N/um because the
% thesis uses base-circle tangential displacement as the torsional DOF.
% Convert to conventional angular stiffness (N m/rad) used by this code.
kThetaSunLin = 10e6;
kThetaRingLin = 57.6e6;
kThetaCarrierLin = 10e6;
p.support.sun.ktheta = kThetaSunLin*p.gear.baseRadiusSun^2;
p.support.ring.ktheta = kThetaRingLin*p.gear.baseRadiusRing^2;
p.support.carrier.ktheta = kThetaCarrierLin*p.gear.carrierTorsionalRadius^2;
p.support.sun.ctheta = torsionalDamping(kThetaSunLin,p.body.sun.J, ...
    p.gear.baseRadiusSun,zetaSupport);
p.support.ring.ctheta = torsionalDamping(kThetaRingLin,p.body.ring.J, ...
    p.gear.baseRadiusRing,zetaSupport);
p.support.carrier.ctheta = torsionalDamping(kThetaCarrierLin,p.body.carrier.J, ...
    p.gear.carrierTorsionalRadius,zetaSupport);

% Mesh damping ratio and transmission-error amplitude are reported.
p.mesh.sunPlanet.zeta = 0.07;
p.mesh.ringPlanet.zeta = 0.07;
p.mesh.sunPlanet.teAmplitude = 20e-6;
p.mesh.ringPlanet.teAmplitude = 20e-6;
% Parker tooth-indexing phase relationship used in Sec. 2.3.2.
p.mesh.sunPlanet.tePhase = mod(p.gear.zs*p.gear.planetPhase0,2*pi);
p.mesh.ringPlanet.tePhase = mod(p.gear.zr*p.gear.planetPhase0,2*pi);

% Standard transverse contact ratios from the reported z, m and 20-deg
% pressure angle. These determine the single/double-tooth overlap pattern.
aSp = p.gear.pitchRadiusSun+p.gear.pitchRadiusPlanet;
aRp = p.gear.pitchRadiusRing-p.gear.pitchRadiusPlanet;
p.tvms.sunPlanet.contactRatio = (...
    sqrt(p.gear.addendumRadiusSun^2-p.gear.baseRadiusSun^2)+ ...
    sqrt(p.gear.addendumRadiusPlanet^2-p.gear.baseRadiusPlanet^2)- ...
    aSp*sin(p.gear.alpha_nominal))/(pi*p.gear.module*cos(p.gear.alpha_nominal));
p.tvms.ringPlanet.contactRatio = (...
    sqrt(p.gear.addendumRadiusPlanet^2-p.gear.baseRadiusPlanet^2)- ...
    sqrt(p.gear.addendumRadiusRing^2-p.gear.baseRadiusRing^2)+ ...
    aRp*sin(p.gear.alpha_nominal))/(pi*p.gear.module*cos(p.gear.alpha_nominal));

% Absolute mean mesh stiffness is calculated in the thesis by the potential
% energy method, but no numeric table/face width is supplied. Keep this as an
% explicit reproduction assumption. The crack/healthy ratio is calculated
% from the thesis straight-root-crack section model.
p.mesh.sunPlanet.kMean = 3.0e8;
% Retained only for the legacy route.  The Peng-Yue route ignores this
% value and builds PR stiffness from Eqs. (2-25)-(2-27).
p.mesh.ringPlanet.kMean = 3.6e8;
p.tvms.profileFloor = 0.18;
p.tvms.profileSin2 = 1.0;
p.tvms.faultTooth = 1;
p.tvms.targetPlanet = 1;
p.tvms.planetRingToothOffset = floor(p.gear.zp/2);
p.tvms.rootCrack.angleDeg = 45;
p.tvms.rootCrack.faceCoverage = 1.0;
p.tvms.rootCrack.lookupPoints = 257;
p.tvms.rootCrack.integrationPoints = 800;
p.tvms.rootCrack.minLigamentFraction = 0.03;
p.tvms.rootCrack.sectionModel = 'pengyue_straight_root';

% Fig. 3-14 study settings.
p.reproduction.fig314.crackDepthMm = 4;
p.reproduction.fig314.crackAnglesDeg = [15 30 45 60 75];
p.reproduction.fig314.maxFrequencyHz = 2000;
p.reproduction.fig314.lowFrequencyMaxHz = 180;
p.reproduction.fig314.reportedHalfBacklashDimensionless = 3;
p.reproduction.fig314.reportedFrictionCoefficient = 0.06;

% The source thesis states b=3 and mu=0.06, but does not report the physical
% displacement scale bc needed to convert dimensionless backlash to metres.
% Friction additionally requires per-contact sliding-arm bookkeeping. The
% first-pass Fig. 3-14 runner therefore keeps these two terms disabled and
% records that limitation instead of inventing hidden values.
p.contact.model = 'linear';
p.contact.friction.enabled = false;
p.contact.friction.mu = 0.06;
p.contact.clearanceSunPlanet = 0;
p.contact.clearanceRingPlanet = 0;

% Fields below document what is exact and what remains inferred.
p.reproduction.source = 'Peng Yue NCEPU MSc thesis, 2025, Fig. 3-14';
p.reproduction.reported = struct( ...
    'sunRpm',800,'outputTorqueNm',50,'teAmplitudeUm',20, ...
    'meshDampingRatio',0.07,'frictionCoefficient',0.06, ...
    'halfBacklashDimensionless',3,'crackDepthMm',4, ...
    'crackAnglesDeg',[15 30 45 60 75]);
p.reproduction.assumption = struct( ...
    'sampleRateHz',p.model.fs, ...
    'faceWidthMm',1e3*p.gear.faceWidthSunPlanet, ...
    'youngsModulusGPa',p.material.gear.E/1e9, ...
    'poissonRatio',p.material.gear.nu, ...
    'meanMeshStiffnessSunPlanetNpm',p.mesh.sunPlanet.kMean, ...
    'meanMeshStiffnessRingPlanetNpm',p.mesh.ringPlanet.kMean, ...
    'balancedInputTorqueNm',p.operating.inputTorque, ...
    'frictionApplied',false,'backlashApplied',false, ...
    'foundationModel','none', ...
    'foundationGeometry',['MISSING: r_int and theta_f for sun/planet; ' ...
        'ring rim/body and fillet geometry'], ...
    'foundationDiagnosticCompliance_m_per_N',NaN, ...
    'constantFoundationDisabled',true, ...
    'prInternalGeometry',['standard-involute angular limits inferred from ' ...
        'reported base/addendum/root radii; requires CAD confirmation']);

% Recompute the working pressure angles and map after all geometry changes.
p.gear.alphaSunPlanet = p.gear.alpha_nominal;
p.gear.alphaRingPlanet = p.gear.alpha_nominal;
p.map = pg_dof_map(p.model.nPlanet);
end

function c = supportDamping(k,m,zeta)
c = 2*zeta*sqrt(k*m);
end

function cTheta = torsionalDamping(kLinear,J,rEquivalent,zeta)
mEquivalent = J/rEquivalent^2;
cLinear = 2*zeta*sqrt(kLinear*mEquivalent);
cTheta = cLinear*rEquivalent^2;
end
