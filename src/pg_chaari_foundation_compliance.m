function [compliance,detail] = pg_chaari_foundation_compliance(member,rContact,p,alphaM,cfg)
%PG_CHAARI_FOUNDATION_COMPLIANCE Chaari/Sainsot fillet-foundation term.
%
% Implements Chaari et al. (2009), Eqs. (3)-(5), with the coefficient
% polynomials of Sainsot et al. (2004).  All lengths are in metres and the
% returned compliance is in m/N.  The inner radius and root half-angle are
% deliberately required inputs: module and tooth count do not determine the
% gear body geometry used by the source paper.
%
% cfg.innerRadius      r_int in Fig. 2 [m]
% cfg.rootHalfAngle    theta_f in Fig. 2 [rad]
% cfg.ufMode           'radial' (default), or explicit 'normal'
% alphaM               working pressure angle alpha_m in Chaari Eq. (3)
%
% For the default radial convention, u_f=max(r_contact-r_f,0), i.e. the
% root-to-contact distance represented in Fig. 2.  This convention and all
% intermediate quantities are returned in detail for auditability.

arguments
    member (1,:) char
    rContact (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
    p struct
    alphaM (1,1) double
    cfg struct
end

if ~(isfinite(alphaM) && alphaM >= 0 && alphaM < pi/2)
    error('pg_chaari_foundation_compliance:PressureAngle', ...
        'alphaM must be finite and satisfy 0 <= alphaM < pi/2.');
end

member = lower(strtrim(member));
switch member
    case 'sun'
        rf = p.gear.rootRadiusSun;
        B = p.gear.faceWidthSunPlanet;
    case 'planet'
        rf = p.gear.rootRadiusPlanet;
        B = p.gear.faceWidthSunPlanet;
    otherwise
        error('pg_chaari_foundation_compliance:Member', ...
            'Unsupported member "%s". Expected sun or planet.',member);
end

E = p.material.gear.E;
if ~(isfinite(E) && E > 0 && isfinite(B) && B > 0 && isfinite(rf) && rf > 0)
    error('pg_chaari_foundation_compliance:MaterialGeometry', ...
        'E, face width, and root radius must be finite positive values.');
end

if ~isfield(cfg,'innerRadius') || ~isfinite(cfg.innerRadius) || cfg.innerRadius <= 0
    error('pg_chaari_foundation_compliance:MissingGeometry', ...
        'Explicit innerRadius (r_int) is required for the %s gear.',member);
end
if ~isfield(cfg,'rootHalfAngle') || ~isfinite(cfg.rootHalfAngle) || cfg.rootHalfAngle <= 0
    error('pg_chaari_foundation_compliance:MissingGeometry', ...
        'Explicit rootHalfAngle (theta_f) is required for the %s gear.',member);
end

rInt = cfg.innerRadius;
thetaF = cfg.rootHalfAngle;
if rInt >= rf
    error('pg_chaari_foundation_compliance:Geometry', ...
        'innerRadius (%.6g m) must be smaller than rootRadius (%.6g m).',rInt,rf);
end
if thetaF >= pi/2
    error('pg_chaari_foundation_compliance:Geometry', ...
        'rootHalfAngle must be in (0,pi/2), got %.6g rad.',thetaF);
end

if isfield(cfg,'ufMode') && ~isempty(cfg.ufMode)
    ufMode = lower(string(cfg.ufMode));
else
    ufMode = "radial";
end
switch ufMode
    case "radial"
        uf = max(rContact-rf,0);
        ufDefinition = 'u_f=max(r_contact-r_f,0) (radial root-to-contact distance)';
    case "normal"
        uf = max(rContact-rf,0)*cos(alphaM);
        ufDefinition = 'u_f=max(r_contact-r_f,0)*cos(alpha_m) (normal projection)';
    otherwise
        error('pg_chaari_foundation_compliance:UfMode', ...
            'Unsupported ufMode "%s". Use "radial" or "normal".',ufMode);
end

% Fig. 2 root geometry: S_f is the chordal thickness at the root circle.
Sf = 2*rf*sin(thetaF);
hfi = rf/rInt;
ratio = uf/Sf;

% Sainsot/Chaari Table 1 coefficients [A B C D E F], rows L*,M*,P*,Q*.
coef = [ ...
    -5.574e-5,  -1.9986e-3, -2.3015e-4,  4.7702e-3,  0.0271, 6.8045; ...
     60.111e-5,  28.100e-3, -83.431e-4, -9.9256e-3,  0.1624,0.9086; ...
    -50.952e-5, 185.50e-3,   0.0538e-4, 53.300e-3,  0.2895,0.9236; ...
     -6.2042e-5,  9.0889e-3, -4.0964e-4,  7.8297e-3, -0.1472,0.6904];
stars = zeros(4,1);
for i = 1:4
    A = coef(i,1); Bi = coef(i,2); C = coef(i,3);
    D = coef(i,4); Ei = coef(i,5); F = coef(i,6);
    stars(i) = A/thetaF^2 + Bi*hfi^2 + C*hfi/thetaF + ...
        D/thetaF + Ei*hfi + F;
end
Lstar = stars(1);
Mstar = stars(2);
Pstar = stars(3);
Qstar = stars(4);

bracket = Lstar*ratio^2 + Mstar*ratio + ...
    Pstar*(1 + Qstar*tan(alphaM)^2);
if ~(isfinite(bracket) && bracket > 0)
    error('pg_chaari_foundation_compliance:NonPositiveCompliance', ...
        ['Chaari polynomial produced a non-positive bracket (%.6g). ' ...
         'Check r_int, theta_f, and u_f geometry.'],bracket);
end

compliance = cos(alphaM)^2/(B*E)*bracket;
if ~(isfinite(compliance) && compliance >= 0)
    error('pg_chaari_foundation_compliance:Numerical', ...
        'Computed foundation compliance is not finite/nonnegative.');
end

detail = struct();
detail.member = member;
detail.rContact = rContact;
detail.rootRadius = rf;
detail.innerRadius = rInt;
detail.rootHalfAngle = thetaF;
detail.hfi = hfi;
detail.uf = uf;
detail.Sf = Sf;
detail.ufOverSf = ratio;
detail.alphaM = alphaM;
detail.faceWidth = B;
detail.E = E;
detail.Lstar = Lstar;
detail.Mstar = Mstar;
detail.Pstar = Pstar;
detail.Qstar = Qstar;
detail.bracket = bracket;
detail.deltaPerUnitForce = compliance;
detail.kFoundation = 1/compliance;
if exist('ufDefinition','var')
    detail.ufDefinition = ufDefinition;
end
detail.formula = ['c_f = cos(alpha_m)^2/(B E) * [' ...
    'L*(u_f/S_f)^2 + M*(u_f/S_f) + P*(1 + Q*tan(alpha_m)^2)]'];
detail.source = 'Chaari et al. 2009 Eqs. (3)-(5); Sainsot et al. 2004 coefficients';
end
