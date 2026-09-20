function [contactFraction,ratio,detail] = pg_pengyue_root_crack_energy_lookup(p,fault)
%PG_PENGYUE_ROOT_CRACK_ENERGY_LOOKUP Crack stiffness ratio from Eqs. (3-1)-(3-3).
% The thesis assumes a straight crack from the tooth root. The crack changes
% the sun-tooth bending and shear sections; Hertz, axial compression and the
% mating planet remain healthy. The absolute pair stiffness continues to
% scale the isolated reproduction's mean-TVMS assumption because tooth face
% width/foundation constants are not fully reported in the thesis.

nLookup = p.tvms.rootCrack.lookupPoints;
contactFraction = linspace(0,1,nLookup).';
ratioFullFace = zeros(nLookup,1);
healthyCompliance = zeros(nLookup,1);
crackedCompliance = zeros(nLookup,1);

for k = 1:nLookup
    fSun = contactFraction(k);
    fPlanet = 1-fSun;
    sunHealthy = toothCompliance('sun',fSun,p,0,0,false);
    sunCracked = toothCompliance('sun',fSun,p, ...
        fault.crackDepth,fault.crackAngle,true);
    planetHealthy = toothCompliance('planet',fPlanet,p,0,0,false);

    kHertz = pi*p.material.gear.E*p.gear.faceWidthSunPlanet/ ...
        (4*(1-p.material.gear.nu^2));
    common = 1/kHertz + sunHealthy.axial + planetHealthy.axial + ...
        planetHealthy.bending + planetHealthy.shear;
    healthyCompliance(k) = common + sunHealthy.bending + sunHealthy.shear;
    crackedCompliance(k) = common + sunCracked.bending + sunCracked.shear;
    ratioFullFace(k) = healthyCompliance(k)/crackedCompliance(k);
end

ratio = min(max(ratioFullFace,p.tvms.rootCrack.minLigamentFraction),1);
detail.fullFaceRatio = ratioFullFace;
detail.healthyCompliance = healthyCompliance;
detail.crackedCompliance = crackedCompliance;
detail.minimumRatio = min(ratio);
detail.maximumDrop = 1-min(ratio);
detail.method = ['Peng Yue straight-root crack: Eqs. (3-1)-(3-3), ' ...
    'sun bending/shear section changed; Hertz/axial retained'];
end

function c = toothCompliance(member,contactFraction,p,crackDepth,crackAngle,isCracked)
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
        error('pg_pengyue_root_crack_energy_lookup:Member','Unknown gear member.');
end

rContact = rb + contactFraction*(ra-rb);
beamLength = max(rContact-rf,1e-6);
nInt = p.tvms.rootCrack.integrationPoints;
dx = beamLength/nInt;
x = ((1:nInt)-0.5)*dx;
r = rf+x;
thicknessHealthy = toothThickness(r,z,rb,rf,p.gear.alpha_nominal);
thickness = thicknessHealthy;

if isCracked && crackDepth > 0
    % Eq. (3-1)/(3-2): h_q = h_c - q*sin(gamma). For sections above h_q,
    % the remaining section is h_c + h_x - q*sin(gamma), where 2*h_x is
    % the healthy tooth thickness and h_c is the root half-thickness.
    hx = thicknessHealthy/2;
    hc = thicknessHealthy(1)/2;
    hq = hc-crackDepth*sin(crackAngle);
    affected = hx>hq;
    crackedSection = hc+hx-crackDepth*sin(crackAngle);
    minSection = p.tvms.rootCrack.minLigamentFraction*thicknessHealthy;
    thickness(affected) = max(crackedSection(affected),minSection(affected));
    thickness = min(thickness,thicknessHealthy);
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
% Thesis keeps axial compression healthy under the crack.
AHealthy = b*thicknessHealthy;
c.axial = sum((sin(alpha)^2./(E*AHealthy))*dx);
end

function thickness = toothThickness(r,z,rb,rf,alpha0)
halfPitchAngle = pi/(2*z);
invAlpha0 = tan(alpha0)-alpha0;
alphaR = zeros(size(r));
aboveBase = r>rb;
alphaR(aboveBase) = acos(rb./r(aboveBase));
halfAngle = halfPitchAngle+invAlpha0-(tan(alphaR)-alphaR);
halfAngle(r<=rb) = halfPitchAngle+invAlpha0;
thickness = 2*r.*halfAngle;
rootThickness = 2*rf*(halfPitchAngle+invAlpha0);
thickness = max(thickness,0.05*rootThickness);
end
