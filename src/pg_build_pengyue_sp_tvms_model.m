function model = pg_build_pengyue_sp_tvms_model(p, fault, options)
%PG_BUILD_PENGYUE_SP_TVMS_MODEL
% Build a precomputed sun-planet (external mesh) TVMS model following the
% potential-energy framework used in Peng Yue's thesis, Eqs. (2-6)-(2-24)
% and the straight root-crack section model, Eqs. (3-1)-(3-3).
%
% IMPORTANT
% 1) This function replaces the old smooth/sin^2 tooth-pair profile.
% 2) Crack changes ONLY sun-tooth bending/shear sections. Hertz contact,
%    axial compression, and the mating planet remain healthy, matching the
%    thesis assumption.
% 3) Tooth-foundation compliance is evaluated with the Chaari/Sainsot
%    because the thesis refers its L*,M*,P*,Q* details to Chaari et al. and
%    does not tabulate all geometry required to recover a unique value.
% 4) A SINGLE global stiffness scale may be supplied for calibration to the
%    healthy Fig. 3-2 curve. Never fit a separate scale for each crack case.
%
% INPUT
%   p     : project parameter struct (pg_parameters_pengyue_fig314)
%   fault : struct with fields
%           type='none' or 'sun'
%           tooth (default 1)
%           crackDepth [m]
%           crackAngle [rad]
%   options.LookupPoints            default 801
%   options.IntegrationPoints       default 1200
%   options.FoundationModel        default "none".  Use "chaari" only
%                                  when the gear-body inner radius and root
%                                  half-angle are supplied explicitly.
%   FoundationCompliance*         retained as compatibility arguments but
%                                  nonzero values are rejected.  The former
%                                  7e-9 m/N constant is diagnostic only and
%                                  is not part of the formal model.
%   options.FoundationInnerRadiusSun/Planet [m].  Required by "chaari"
%                                  unless supplied in p.gear.foundation.*.
%   options.FoundationRootHalfAngleSun/Planet [rad].  Required by "chaari"
%                                  unless supplied in p.gear.foundation.*.
%   options.FoundationUfMode       default "radial".  The Chaari Fig. 2
%                                  definition is represented by the radial
%                                  root-to-contact distance r_c-r_f; this is
%                                  explicit so it cannot be mistaken for a
%                                  fitted scale.
%   options.GlobalScale             default 1
%
% OUTPUT model
%   .eta                 pair contact progress 0..1
%   .kPairHealthy        healthy single-pair stiffness [N/m]
%   .kPairFault          cracked single-pair stiffness [N/m]
%   .contactRatio        transverse contact ratio
%   .meshPhaseCycles     Parker phase in mesh-cycle units for each planet
%   .globalScale         one common absolute stiffness scale
%
% Usage:
%   p = pg_parameters_pengyue_fig314;
%   f = struct('type','sun','tooth',1,'crackDepth',3e-3,...
%              'crackAngle',deg2rad(45));
%   mdl = pg_build_pengyue_sp_tvms_model(p,f);
%   [k,detail] = pg_eval_pengyue_sp_tvms(t,p,f,mdl);

arguments
    p struct
    fault struct
    options.LookupPoints (1,1) double {mustBeInteger,mustBeGreaterThan(options.LookupPoints,100)} = 801
    options.IntegrationPoints (1,1) double {mustBeInteger,mustBeGreaterThan(options.IntegrationPoints,200)} = 1200
    options.FoundationModel (1,1) string = "none"
    options.FoundationComplianceSun (1,1) double {mustBeNonnegative} = 0
    options.FoundationCompliancePlanet (1,1) double {mustBeNonnegative} = 0
    options.FoundationInnerRadiusSun (1,1) double = NaN
    options.FoundationInnerRadiusPlanet (1,1) double = NaN
    options.FoundationRootHalfAngleSun (1,1) double = NaN
    options.FoundationRootHalfAnglePlanet (1,1) double = NaN
    options.FoundationUfMode (1,1) string = "missing"
    options.GlobalScale (1,1) double {mustBePositive} = 1
end

foundationModel = lower(strtrim(options.FoundationModel));
if abs(options.GlobalScale-1) > 10*eps
    error('pg_build_pengyue_sp_tvms_model:GlobalScaleForbidden', ...
        'GlobalScale is disabled for the physical reproduction route; use 1.');
end
if ~ismember(foundationModel,["none","chaari"])
    error('pg_build_pengyue_sp_tvms_model:FoundationModel', ...
        'FoundationModel must be "none" or "chaari".');
end
if options.FoundationComplianceSun > 0 || options.FoundationCompliancePlanet > 0
    error('pg_build_pengyue_sp_tvms_model:ConstantFoundationForbidden', ...
        ['Nonzero FoundationCompliance* belongs to the retired diagnostic ' ...
         'route. Supply explicit Chaari geometry or use zero with FoundationModel="none".']);
end
if foundationModel == "chaari"
    foundationCfg.sun = resolveFoundationGeometry(p,'sun', ...
        options.FoundationInnerRadiusSun,options.FoundationRootHalfAngleSun, ...
        options.FoundationUfMode);
    foundationCfg.planet = resolveFoundationGeometry(p,'planet', ...
        options.FoundationInnerRadiusPlanet,options.FoundationRootHalfAnglePlanet, ...
        options.FoundationUfMode);
else
    foundationCfg = struct();
end

% ----- geometry reported / implied by the thesis -----
zs = p.gear.zs;
zp = p.gear.zp;
m  = p.gear.module;
alpha0 = p.gear.alpha_nominal;

rs = m*zs/2;
rp = m*zp/2;
rbS = p.gear.baseRadiusSun;
rbP = p.gear.baseRadiusPlanet;
raS = p.gear.addendumRadiusSun;
raP = p.gear.addendumRadiusPlanet;
rfS = p.gear.rootRadiusSun;
rfP = p.gear.rootRadiusPlanet;

% Standard spur-gear transverse contact ratio from the line of action.
pathApproach = sqrt(max(raP^2-rbP^2,0)) - rp*sin(alpha0);
pathRecess   = sqrt(max(raS^2-rbS^2,0)) - rs*sin(alpha0);
basePitch    = pi*m*cos(alpha0);
pathTotal    = pathApproach + pathRecess;
epsilon      = pathTotal/basePitch;

if ~(epsilon > 1 && epsilon < 2)
    warning('pg_build_pengyue_sp_tvms_model:ContactRatio', ...
        'Computed contact ratio %.4f is outside the expected 1..2 range.',epsilon);
end

eta = linspace(0,1,options.LookupPoints).';
kHealthy = zeros(size(eta));
kFault   = zeros(size(eta));
% Columns are [Hertz, sun-bending, sun-axial, sun-shear, sun-foundation,
% planet-bending, planet-axial, planet-shear, planet-foundation, total].
compHealthy = zeros(numel(eta),10);
compFault   = zeros(numel(eta),10);
foundationDetailSun = repmat(emptyFoundationDetail(),numel(eta),1);
foundationDetailPlanet = repmat(emptyFoundationDetail(),numel(eta),1);

for i = 1:numel(eta)
    % Coordinate along line of action, zero at pitch point.
    s = -pathApproach + eta(i)*pathTotal;

    % Distance from base-circle tangent point to the contact point.
    lS = rs*sin(alpha0) + s;
    lP = rp*sin(alpha0) - s;
    lS = max(lS,0);
    lP = max(lP,0);

    rFS = sqrt(rbS^2 + lS^2);
    rFP = sqrt(rbP^2 + lP^2);
    rFS = min(max(rFS,rbS),raS);
    rFP = min(max(rFP,rbP),raP);

    sunH = toothComplianceExternal('sun',rFS,p,false,fault,options.IntegrationPoints);
    sunC = toothComplianceExternal('sun',rFS,p,true, fault,options.IntegrationPoints);
    plaH = toothComplianceExternal('planet',rFP,p,false,fault,options.IntegrationPoints);

    % Eq. (2-21): Hertz contact stiffness.
    E = p.material.gear.E;
    nu = p.material.gear.nu;
    B = p.gear.faceWidthSunPlanet;
    kH = pi*E*B/(4*(1-nu^2));
    cH = 1/kH;

    % Chaari/Sainsot fillet-foundation compliance.  This is deliberately
    % evaluated at the current contact position; it is not a global scale
    % factor and it does not alter contact-ratio or active-pair logic.
    if foundationModel == "chaari"
        alphaM = p.gear.alpha_nominal;
        if isfield(p.gear,'alphaSunPlanet') && isfinite(p.gear.alphaSunPlanet)
            alphaM = p.gear.alphaSunPlanet;
        end
        [cFS,foundationDetailSun(i)] = pg_chaari_foundation_compliance( ...
            'sun',rFS,p,alphaM,foundationCfg.sun);
        [cFP,foundationDetailPlanet(i)] = pg_chaari_foundation_compliance( ...
            'planet',rFP,p,alphaM,foundationCfg.planet);
    else
        cFS = 0;
        cFP = 0;
    end

    % Eq. (2-23), written as total compliance of one tooth pair.
    cHealthy = cH + ...
        sunH.bending + sunH.axial + sunH.shear + cFS + ...
        plaH.bending + plaH.axial + plaH.shear + cFP;

    % Peng Yue Eqs. (3-1)-(3-3): root crack alters sun bending/shear only.
    % Hertz, axial compression and the mating planet remain healthy.
    cFault = cH + ...
        sunC.bending + sunH.axial + sunC.shear + cFS + ...
        plaH.bending + plaH.axial + plaH.shear + cFP;

    kHealthy(i) = 1/cHealthy;
    kFault(i)   = 1/cFault;

    compHealthy(i,:) = [cH sunH.bending sunH.axial sunH.shear cFS ...
        plaH.bending plaH.axial plaH.shear cFP cHealthy];
    compFault(i,:) = [cH sunC.bending sunH.axial sunC.shear cFS ...
        plaH.bending plaH.axial plaH.shear cFP cFault];
end

kHealthy = options.GlobalScale*kHealthy;
kFault   = options.GlobalScale*kFault;

% Parker mesh phase relation gamma_sp^i = z_s * phi_p^i. Convert rad to
% mesh-cycle offsets so each planet uses its own local mesh counter.
planetPhase = p.gear.planetPhase0(:).';
meshPhaseCycles = mod(zs*planetPhase/(2*pi),1);
% For tooth identity, the full (unwrapped) phase is needed as well.
toothPhaseCycles = zs*planetPhase/(2*pi);

model.eta = eta;
model.kPairHealthy = kHealthy;
model.kPairFault = kFault;
model.contactRatio = epsilon;
model.pathApproach = pathApproach;
model.pathRecess = pathRecess;
model.pathTotal = pathTotal;
model.basePitch = basePitch;
model.meshPhaseCycles = meshPhaseCycles;
model.toothPhaseCycles = toothPhaseCycles;
model.globalScale = options.GlobalScale;
model.foundationComplianceSun = options.FoundationComplianceSun;
model.foundationCompliancePlanet = options.FoundationCompliancePlanet;
model.foundationModel = foundationModel;
model.foundationConfig = foundationCfg;
model.foundationDetailSun = foundationDetailSun;
model.foundationDetailPlanet = foundationDetailPlanet;
model.complianceColumnNames = {'hertz','sunBending','sunAxial','sunShear', ...
    'sunFoundation','planetBending','planetAxial','planetShear', ...
    'planetFoundation','total'};
model.complianceHealthy = compHealthy;
model.complianceFault = compFault;
model.method = ['Potential-energy external-mesh TVMS. Numerical tooth-section ' ...
    'integration implements Eqs. (2-6)-(2-10); crack section follows ' ...
    'Eqs. (3-1)-(3-3). Foundation term uses Chaari/Sainsot Eq. (3)-(5) ' ...
    'when FoundationModel="chaari". Single/double pairs are summed in ' ...
    'pg_eval_pengyue_sp_tvms.'];
end

function c = toothComplianceExternal(member,rF,p,applyCrack,fault,nInt)
% Numerical form of Eqs. (2-6)-(2-8). Geometry is represented by an
% involute tooth section. This avoids the old empirical sin^2 profile while
% retaining the same physical strain-energy terms as Peng Yue.

switch lower(member)
    case 'sun'
        z  = p.gear.zs;
        rb = p.gear.baseRadiusSun;
        rf = p.gear.rootRadiusSun;
    case 'planet'
        z  = p.gear.zp;
        rb = p.gear.baseRadiusPlanet;
        rf = p.gear.rootRadiusPlanet;
    otherwise
        error('Unknown member: %s',member);
end

alpha0 = p.gear.alpha_nominal;
E = p.material.gear.E;
G = p.material.gear.G;
B = p.gear.faceWidthSunPlanet;

% Base-circle half tooth angle phi_b used in Eq. (2-17).
phiB = pi/(2*z) + invInvolute(alpha0);

% Use radial stations, then integrate along the tooth-centre-line coordinate.
r = linspace(rf,rF,nInt).';
alphaR = zeros(size(r));
mask = r > rb;
alphaR(mask) = acos(rb./r(mask));

% Involute half-tooth angle at each radius. Below the base circle the model
% continues the base/root section with constant angular half-thickness,
% consistent with the piecewise form in Eq. (2-17).
thetaHalf = phiB - invInvolute(alphaR);
thetaHalf(~mask) = phiB;

xAxis = r.*cos(thetaHalf);
xRoot = rf*cos(phiB);
x = xAxis-xRoot;
hx = r.*sin(thetaHalf);      % healthy half thickness
sHealthy = 2*hx;

% Guard only numerical degeneracy; do not smooth the TVMS waveform.
minHealthy = max(1e-9,1e-4*max(sHealthy));
sHealthy = max(sHealthy,minHealthy);

sBS = sHealthy; % section used by bending/shear
if applyCrack && isfield(fault,'type') && strcmpi(fault.type,'sun') && ...
        isfield(fault,'crackDepth') && fault.crackDepth > 0
    q = fault.crackDepth;
    gamma = fault.crackAngle;

    % Eqs. (3-1),(3-2): h_q = h_c-q*sin(gamma), and for h_x>h_q the
    % remaining full section is h_x+h_c-q*sin(gamma).
    hc = rf*sin(phiB);
    hq = hc-q*sin(gamma);
    affected = hx > hq;
    sCrack = hx + hc - q*sin(gamma);

    % Numerical protection only. Keep the ligament limit very small so that
    % it does not act as an artificial stiffness-shaping function.
    minLig = 0.005*sHealthy;
    sBS(affected) = max(sCrack(affected),minLig(affected));
    sBS = min(sBS,sHealthy);
end

Ahealthy = B*sHealthy;
ABS = B*sBS;
IBS = B*sBS.^3/12;

alphaF = acos(min(max(rb/rF,0),1));
xF = x(end);
hF = hx(end);
FaOverF = sin(alphaF);  % Eq. (2-13)
FbOverF = cos(alphaF);  % Eq. (2-14)

% Eq. (2-6): M/F = Fb/F*(d-x)-Fa/F*h.
MOverF = FbOverF.*(xF-x) - FaOverF.*hF;

% Remove duplicate/non-monotonic coordinates before trapz.
[xu,ia] = unique(x,'stable');
MOverF = MOverF(ia);
IBS = IBS(ia);
ABS = ABS(ia);
Ahealthy = Ahealthy(ia);

c.bending = trapz(xu, MOverF.^2./(E*IBS));
% Peng Yue states axial compression remains healthy for the crack case.
c.axial = trapz(xu, (FaOverF^2)./(E*Ahealthy));
c.shear = trapz(xu, (1.2*FbOverF^2)./(G*ABS));
c.alphaF = alphaF;
c.hF = hF;
c.xF = xF;
end

function y = invInvolute(alpha)
y = tan(alpha)-alpha;
end

function cfg = resolveFoundationGeometry(p,member,innerRadius,rootHalfAngle,ufMode)
% Resolve only explicitly supplied Chaari geometry.  No inner radius or
% root-angle value is inferred from module/tooth count here.
cfg = struct('innerRadius',innerRadius,'rootHalfAngle',rootHalfAngle, ...
    'ufMode',ufMode);
if lower(strtrim(ufMode)) == "missing"
    error('pg_build_pengyue_sp_tvms_model:MissingFoundationUfDefinition', ...
        ['FoundationModel="chaari" also requires an explicit u_f geometry ' ...
         'definition. Select a documented mode only after confirming the ' ...
         'root/contact coordinates.']);
end
if ~isfinite(cfg.innerRadius)
    if isfield(p,'gear') && isfield(p.gear,'foundation') && ...
            isfield(p.gear.foundation,member) && ...
            isfield(p.gear.foundation.(member),'innerRadius')
        cfg.innerRadius = p.gear.foundation.(member).innerRadius;
    end
end
if ~isfinite(cfg.rootHalfAngle)
    if isfield(p,'gear') && isfield(p.gear,'foundation') && ...
            isfield(p.gear.foundation,member) && ...
            isfield(p.gear.foundation.(member),'rootHalfAngle')
        cfg.rootHalfAngle = p.gear.foundation.(member).rootHalfAngle;
    end
end
if ~isfinite(cfg.innerRadius) || ~isfinite(cfg.rootHalfAngle)
    missing = {};
    if ~isfinite(cfg.innerRadius); missing{end+1} = 'innerRadius'; end %#ok<AGROW>
    if ~isfinite(cfg.rootHalfAngle); missing{end+1} = 'rootHalfAngle'; end %#ok<AGROW>
    error('pg_build_pengyue_sp_tvms_model:MissingFoundationGeometry', ...
        ['FoundationModel="chaari" requires explicit %s geometry for ' ...
         '%s gear. No value is inferred from module/tooth count.'], ...
        strjoin(missing,', '),member);
end
end

function d = emptyFoundationDetail()
% Stable template for model.foundationDetail* when the model is disabled.
d = struct('member','','rContact',NaN,'rootRadius',NaN,'innerRadius',NaN, ...
    'rootHalfAngle',NaN,'hfi',NaN,'uf',NaN,'Sf',NaN,'ufOverSf',NaN, ...
    'alphaM',NaN,'faceWidth',NaN,'E',NaN,'Lstar',NaN,'Mstar',NaN, ...
    'Pstar',NaN,'Qstar',NaN,'bracket',NaN,'deltaPerUnitForce',NaN, ...
    'kFoundation',NaN,'ufDefinition','','formula','','source','');
end
