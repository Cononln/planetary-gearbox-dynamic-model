function p = pg_parameters
%PG_PARAMETERS Frozen geometry and initial calibration parameters in SI.
% Values labelled estimated must later be identified or sensitivity-tested.

p.model.name = 'ZLS160 three-planet 18-DOF baseline';
p.model.nPlanet = 3;
p.model.fs = 51200;
p.model.phi_c0 = 0;

p.gear.zs = 21;
p.gear.zp = 31;
p.gear.zr = 84;
p.gear.module = 1.5e-3;
p.gear.alpha_nominal = deg2rad(20);
% User-confirmed effective face widths for both mesh pairs.
p.gear.faceWidthSunPlanet = 21.5e-3;
p.gear.faceWidthRing = 21.5e-3;
% User-selected nominal working centre for the current model.  It is a
% modelling convention pending direct assembly measurement, not a CAD fact.
p.gear.carrierRadius = 39.75e-3;

p.gear.pitchRadiusSun = p.gear.module*p.gear.zs/2;
p.gear.pitchRadiusPlanet = p.gear.module*p.gear.zp/2;
p.gear.pitchRadiusRing = p.gear.module*p.gear.zr/2;
p.gear.baseRadiusSun = p.gear.pitchRadiusSun*cos(p.gear.alpha_nominal);
p.gear.baseRadiusPlanet = p.gear.pitchRadiusPlanet*cos(p.gear.alpha_nominal);
p.gear.baseRadiusRing = p.gear.pitchRadiusRing*cos(p.gear.alpha_nominal);
% User-confirmed physical outside/root dimensions of the external gears;
% reported inner tooth-tip/root radii of the fixed ring.  They are kept
% separate from pitch/base geometry: dimensions alone do not identify the
% manufacturing profile shifts, tooth thickness, or root fillet.
p.gear.addendumRadiusSun = 34.50e-3/2;
p.gear.rootRadiusSun = 27.75e-3/2;
p.gear.addendumRadiusPlanet = 49.50e-3/2;
p.gear.rootRadiusPlanet = 42.75e-3/2;
p.gear.addendumRadiusRing = 61.50e-3; % internal tooth-tip radius
% Ring dimensions extracted from the supplied CAD geometry (2026-09-08).
p.gear.rootRadiusRing = 64.875e-3;   % internal tooth-root radius
p.gear.outerRadiusRing = 81.000e-3; % housing-integral outer rim radius
p.gear.ringHoleRadius = 72.500e-3;  % eight through-hole pitch radius
p.gear.ringHoleDiameter = 6.8e-3;
p.gear.ringHoleCount = 8;

% User-specified generic structural steel.  These are engineering defaults,
% not a material-certificate identification of the physical gears.
p.material.gear.E = 200e9;
p.material.gear.nu = 0.30;
p.material.gear.rho = 7850;
p.material.gear.G = p.material.gear.E/(2*(1+p.material.gear.nu));

a0sp = p.gear.pitchRadiusSun + p.gear.pitchRadiusPlanet;
cosAlphaSp = a0sp*cos(p.gear.alpha_nominal)/p.gear.carrierRadius;
p.gear.alphaSunPlanet = acos(cosAlphaSp);
p.gear.alphaRingPlanet = p.gear.alpha_nominal;
p.gear.planetPhase0 = 2*pi*(0:p.model.nPlanet-1)/p.model.nPlanet;

% CAD mass and inertia about the rotating z axis.
p.body.sun.mass = 2.00004157544214;
p.body.sun.J = 8.18267837505362e-4;
p.body.planet.mass = 0.206246750168060;
p.body.planet.J = 7.11954055833396e-5;
p.body.carrier.mass = 3.2889145;
p.body.carrier.J = 2.714205e-3;
p.body.ring.mass = 3.21735309060989;
p.body.ring.J = 1.70127867386670e-2;

p.operating.sunRpm = 600;
p.operating.inputTorque = 0;
p.operating.outputTorque = 0;
p.operating.torqueRippleFraction = 0;
p.operating.torqueRippleFrequencyHz = 10;
p.operating.encoderPpr = 1024;
p.kin = pg_kinematics(p.gear.zs, p.gear.zp, p.gear.zr, ...
                      p.operating.sunRpm, p.operating.encoderPpr);
p.operating.torqueRippleFrequencyHz = p.kin.f_s;

% Initial bearing/support values: calibration parameters, not CAD facts.
p.support.sun.kxy = 1.2e8;       p.support.sun.cxy = 9.0e2;
p.support.sun.ktheta = 3.0e3;    p.support.sun.ctheta = 1.0;
p.support.ring.kxy = 3.0e8;      p.support.ring.cxy = 1.8e3;
p.support.ring.ktheta = 2.0e6;   p.support.ring.ctheta = 1.5e2;
p.support.carrier.kxy = 8.0e7;   p.support.carrier.cxy = 8.0e2;
p.support.carrier.ktheta = 0;    p.support.carrier.ctheta = 0.5;
p.support.planetPin.kxy = 1.5e8; p.support.planetPin.cxy = 7.0e2;

% Stage-1 constant mean mesh stiffness. Stage 2 replaces these by TVMS.
p.mesh.sunPlanet.kMean = 4.5e8;
p.mesh.ringPlanet.kMean = 6.0e8;
p.mesh.sunPlanet.zeta = 0.03;
p.mesh.ringPlanet.zeta = 0.03;

% Small deterministic transmission errors for the baseline forcing interface.
p.mesh.sunPlanet.teAmplitude = 0.5e-6;
p.mesh.ringPlanet.teAmplitude = 0.5e-6;
% Mesh phase is set by tooth indexing at each equally spaced planet, not by
% the planets' geometric 120-deg separation itself.  Here zs=21 and zr=84
% are both divisible by three, so all three meshes are in phase.  Computing
% the values from the tooth counts keeps that condition explicit and avoids
% the artificial 166/172-Hz lines produced by [0,120,240] deg TE phases.
p.mesh.sunPlanet.tePhase = mod(p.gear.zs*p.gear.planetPhase0,2*pi);
p.mesh.ringPlanet.tePhase = mod(p.gear.zr*p.gear.planetPhase0,2*pi);
p.mesh.sunPlanet.tePhase(abs(p.mesh.sunPlanet.tePhase)<1e-12) = 0;
p.mesh.ringPlanet.tePhase(abs(p.mesh.ringPlanet.tePhase)<1e-12) = 0;

% Contact ratios computed from the confirmed gear geometry.
p.tvms.sunPlanet.contactRatio = 1.151;
p.tvms.ringPlanet.contactRatio = 1.934;
p.tvms.nWidthSlices = 80;
% Smooth tooth-pair load sharing: zero contribution at contact entry/exit.
% The total mesh stiffness remains nonzero because contactRatio > 1.
p.tvms.profileFloor = 0;
p.tvms.profileSin2 = 1;
p.tvms.faultTooth = 1;
p.tvms.targetPlanet = 1;
p.tvms.planetRingToothOffset = floor(p.gear.zp/2);
p.tvms.radialRemovalDepth = 3.375e-3;
% Root-crack defaults follow Shen (2023). Crack depth remains an explicit
% case parameter; no unreported speed-ripple amplitude is assumed.
p.tvms.rootCrack.angleDeg = 65;
p.tvms.rootCrack.faceCoverage = 1.0;
p.tvms.rootCrack.lookupPoints = 257;
p.tvms.rootCrack.integrationPoints = 600;
p.tvms.rootCrack.minLigamentFraction = 0.05;

% Analytical ring-mode damping and radial force participation are calibrated
% quantities. Frequencies themselves are loaded from the supplied table.
p.ringModal.zeta = 0.010;
p.ringModal.waveNumbers = [2,2,3,3,4,4];
p.ringModal.family = {'cos','sin','cos','sin','cos','sin'};
p.ringModal.sourceModes = 7:12;
p.ringModal.radialMeshProjection = -sin(p.gear.alphaRingPlanet);

% Default sensor model for the pure lumped-parameter workflow.  The
% 18-DOF gearbox generates the actual mesh-contact forces.  Each fixed
% accelerometer is represented by a one-way local mass-spring-damper
% observer driven through a smooth moving-source proximity window.  This
% retains the time-varying planet-to-sensor path without an FE/modal ring.
% Sensor geometry is split into installation position and sensitive-axis
% direction. They are identical for the present two radial accelerometers,
% but keeping them separate supports oblique sensors without changing the
% gearbox equations. `angles` remains a compatibility alias.
p.sensor.locationAngles = [0,pi/2];
p.sensor.sensitiveAngles = [0,pi/2];
p.sensor.angles = p.sensor.locationAngles;
p.sensor.pathFloor = 0.08;
p.sensor.pathKappa = 2.5;
% Small zero-mean source/path asymmetry represents unequal planet loading,
% pin-position tolerances, and non-axisymmetric casing transmission.  With
% [1,1,1] the fc and 2*fc components cancel exactly; the present +/-4%
% engineering estimate restores the carrier-order sidebands observed in
% real gearboxes.  Replace these factors after load-sharing identification.
p.sensor.planetSourceGain = [1.00,0.96,1.04];
p.sensor.pathEffectiveMass = p.body.ring.mass/2;
p.sensor.pathNaturalFrequencyHz = mean([1832.565469,1845.781068]);
p.sensor.pathZeta = 0.02;
p.sensor.rigidRingContribution = 1.0;

% Optional nonlinear tooth-contact law. Validated paper figures retain the
% linear model by default. Use 'unilateral' or 'backlash' only for a
% calibrated contact-loss/clearance sensitivity study.
p.contact.model = 'linear';
p.contact.clearanceSunPlanet = 0;
p.contact.clearanceRingPlanet = 0;
p.contact.dampingOnlyInContact = true;
p.contact.friction.enabled = false;
p.contact.friction.mu = 0.06;

p.map = pg_dof_map(p.model.nPlanet);
end
