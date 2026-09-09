function c = tvms_config(p)
%TVMS_CONFIG Inputs for healthy potential-energy TVMS.
% Values copied from pg_parameters are source values.  NaN entries are
% deliberately unresolved: do not replace them with standard-gear defaults.

c.source = 'src/pg_parameters.m + docs/CURRENT_STATE.md';
c.nPhase = 2001;                 % >= 2000 points over [0,1]
c.outputUnit = 'N/m';
c.allowUnresolvedGeometry = false;

c.zs = p.gear.zs; c.zp = p.gear.zp; c.zr = p.gear.zr;
c.module = p.gear.module;
c.alphaNominal = p.gear.alpha_nominal;
c.alphaSP = p.gear.alphaSunPlanet;
c.alphaPR = p.gear.alphaRingPlanet;
c.centerSP = p.gear.carrierRadius;
c.centerPR = p.gear.carrierRadius;
c.faceWidthSP = p.gear.faceWidthSunPlanet;
c.faceWidthPR = p.gear.faceWidthRing;
c.E = p.material.gear.E; c.nu = p.material.gear.nu;
c.rho = p.material.gear.rho;

% Measured external-gear dimensions present in pg_parameters.m.
c.sun.tipRadius = p.gear.addendumRadiusSun;
c.sun.rootRadius = p.gear.rootRadiusSun;
c.sun.addendumCoefficient = 1;   % h in the supplied sun drawing
c.sun.commonNormalLength = 11.512e-3; % drawing: W across k=3 teeth
c.sun.commonNormalSpanTeeth = 3;
c.planet.tipRadius = p.gear.addendumRadiusPlanet;
c.planet.rootRadius = p.gear.rootRadiusPlanet;
c.planet.addendumCoefficient = 1; % h in the supplied planet drawing
c.planet.commonNormalLength = 16.15e-3; % drawing: W across k=4 teeth
c.planet.commonNormalSpanTeeth = 4;

% Reported physical internal-tooth dimensions.  They define the contact
% boundary but do not identify the root fillet or the profile shift.
c.ring.tipRadius = p.gear.addendumRadiusRing;
c.ring.rootRadius = p.gear.rootRadiusRing;
c.haStar = 1;
c.clearanceStar = 0.25;
% The supplied sun and planet drawings explicitly state xn=0.  For beta=0,
% this is the transverse profile shift used here.  Ring shift remains open.
c.xSun = 0;
c.xPlanet = 0;
c.xRing = 0;

% User-authorised modelling assumption (2026-09-07): retain the present
% 39.75 mm SP working centre and use the implied correction to tooth
% thickness.  It is not claimed to be the manufactured profile shift.
c.useEquivalentWorkingCentreThicknessCorrection = true;
c.equivalentCorrectionAllocation = 'symmetric_sun_planet';

% Root/foundation compliance needs an identified model or measured root
% geometry.  This is intentionally not an empirical tuneable constant.
c.foundationModel = 'annular_rim_plane_stress_fem';
% User-provided ring geometry: 16.125 mm nominal rim thickness, eight
% through holes and an outer rim integral with the gearbox housing.  The
% Python evaluator uses a 2-D plane-stress annulus with the outer boundary
% fixed; the unknown tooth fillet remains an explicit approximation.
c.ring.outerRadius = p.gear.outerRadiusRing;
c.ring.holeRadius = p.gear.ringHoleRadius;
c.ring.holeDiameter = p.gear.ringHoleDiameter;
c.ring.holeCount = p.gear.ringHoleCount;
c.ring.outerBoundary = 'housing_integral_clamped';
% No exact root-fillet curve is available; the root load patch uses the
% involute tooth thickness at the reported root radius.
c.foundationCalibration = 1;
end
