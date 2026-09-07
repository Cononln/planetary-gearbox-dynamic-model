function [contactFraction,ratio,detail] = pg_root_crack_energy_lookup(p,fault)
%PG_ROOT_CRACK_ENERGY_LOOKUP Relative mesh stiffness for a sun-root crack.
% The healthy and cracked tooth pairs share unchanged Hertz and axial
% compression terms. Only the sun-tooth bending and shear compliances use
% the crack-reduced section A(x), I(x), matching the modeling assumption in
% Shen (2023), Eqs. (2.1)-(2.15). The analytical ratio scales the existing
% calibrated tooth-pair stiffness rather than replacing its absolute level.

nLookup = p.tvms.rootCrack.lookupPoints;
contactFraction = linspace(0,1,nLookup).';
ratioFullFace = zeros(nLookup,1);
healthyCompliance = zeros(nLookup,1);
crackedCompliance = zeros(nLookup,1);

for k = 1:nLookup
    fSun = contactFraction(k);
    fPlanet = 1-fSun;
    sunHealthy = toothCompliance('sun',fSun,p,0,0);
    sunCracked = toothCompliance('sun',fSun,p, ...
        fault.crackDepth,fault.crackAngle);
    planetHealthy = toothCompliance('planet',fPlanet,p,0,0);

    % Hertz and axial compression remain unaffected by the root crack.
    kHertz = pi*p.material.gear.E*p.gear.faceWidthSunPlanet/ ...
        (4*(1-p.material.gear.nu^2));
    common = 1/kHertz + sunHealthy.axial + planetHealthy.axial + ...
        planetHealthy.bending + planetHealthy.shear;
    healthyCompliance(k) = common + sunHealthy.bending + sunHealthy.shear;
    crackedCompliance(k) = common + sunCracked.bending + sunCracked.shear;
    ratioFullFace(k) = healthyCompliance(k)/crackedCompliance(k);
end

% Across the face width, intact and cracked strips act in parallel. This
% also permits the measured crack span to be used later without changing
% the two-dimensional energy model.
ratio = (1-fault.faceCoverage) + fault.faceCoverage*ratioFullFace;
ratio = min(max(ratio,p.tvms.rootCrack.minLigamentFraction),1);

detail.fullFaceRatio = ratioFullFace;
detail.healthyCompliance = healthyCompliance;
detail.crackedCompliance = crackedCompliance;
detail.minimumRatio = min(ratio);
detail.maximumDrop = 1-min(ratio);
detail.method = ['Energy-equivalent tapered involute tooth; crack reduces ' ...
    'sun bending/shear section only; Hertz/axial terms unchanged'];
end

function c = toothCompliance(member,contactFraction,p,crackDepth,crackAngle)
switch member
    case 'sun'
        z = p.gear.zs;
        rb = p.gear.baseRadiusSun;
        rf = p.gear.rootRadiusSun;
        ra = p.gear.addendumRadiusSun;
    case 'planet'
        z = p.gear.zp;
        rb = p.gear.baseRadiusPlanet;
        rf = p.gear.rootRadiusPlanet;
        ra = p.gear.addendumRadiusPlanet;
    otherwise
        error('pg_root_crack_energy_lookup:Member','Unknown gear member.');
end

% Contact progresses from the sun base-side to addendum-side. The planet
% fraction is reversed by the caller. Root-to-base material is retained in
% the cantilever integral even though the involute begins at rb.
rContact = rb + contactFraction*(ra-rb);
beamLength = max(rContact-rf,1e-6);
nInt = p.tvms.rootCrack.integrationPoints;
dx = beamLength/nInt;
x = ((1:nInt)-0.5)*dx;
r = rf+x;
thicknessHealthy = toothThickness(r,z,rb,rf,p.gear.alpha_nominal);
thickness = thicknessHealthy;

if crackDepth > 0
    normalPenetration = crackDepth*sin(crackAngle);
    radialReach = max(crackDepth*cos(crackAngle),dx);
    affected = x <= radialReach;
    minThickness = p.tvms.rootCrack.minLigamentFraction* ...
        thicknessHealthy(affected);
    thickness(affected) = max(thicknessHealthy(affected)-normalPenetration, ...
        minThickness);
end

b = p.gear.faceWidthSunPlanet;
A = b*thickness;
I = b*thickness.^3/12;
alpha = p.gear.alpha_nominal;
lever = (beamLength-x)*cos(alpha);
E = p.material.gear.E;
G = p.material.gear.G;

c.bending = sum((lever.^2./(E*I))*dx);
c.shear = sum((1.2*cos(alpha)^2./(G*A))*dx);
% Retained from the healthy section for both healthy and cracked calls.
AHealthy = b*thicknessHealthy;
c.axial = sum((sin(alpha)^2./(E*AHealthy))*dx);
end

function thickness = toothThickness(r,z,rb,rf,alpha0)
halfPitchAngle = pi/(2*z);
invAlpha0 = tan(alpha0)-alpha0;
alphaR = zeros(size(r));
aboveBase = r>rb;
alphaR(aboveBase) = acos(rb./r(aboveBase));
halfAngle = halfPitchAngle + invAlpha0 - (tan(alphaR)-alphaR);
halfAngle(r<=rb) = halfPitchAngle+invAlpha0;
thickness = 2*r.*halfAngle;

% The root fillet is not supplied as a curve. Continue the base-circle
% tooth angle to rf, which is explicit and reproducible, and guard only
% against numerical degeneracy near the addendum.
rootThickness = 2*rf*(halfPitchAngle+invAlpha0);
thickness = max(thickness,0.05*rootThickness);
end
