function model = pg_build_pengyue_pr_tvms_model(p, fault, options)
%PG_BUILD_PENGYUE_PR_TVMS_MODEL Physical planet--ring TVMS lookup.
%
% This is the internal-mesh counterpart of
% pg_build_pengyue_sp_tvms_model.  The planet tooth uses the same numerical
% potential-energy section integration as the validated SP model.  The ring
% tooth uses Peng Yue Eqs. (2-25)--(2-27), i.e. the internal-gear bending,
% radial-compression and shear integrals.  Active tooth pairs are summed in
% parallel by pg_eval_pengyue_pr_tvms.
%
% IMPORTANT GEOMETRY STATUS
% The thesis does not report the ring rim/body geometry required for a
% Chaari/Sainsot ring foundation term.  Therefore the formal default is
% FoundationModel="none" and the ring foundation compliance is exactly
% zero (kf=Inf), with an explicit MISSING status in the model metadata.
% The old 7e-9 m/N diagnostic constant is deliberately not accepted.
%
% The internal-involute formula requires angular limits that are not
% tabulated in the thesis.  When they are absent, this function derives a
% reproducible standard-involute approximation from the reported base,
% addendum and root radii, and labels it as an assumption.  Supplying
% p.gear.internalMeshGeometry overrides those derived values.

arguments
    p struct
    fault struct
    options.LookupPoints (1,1) double {mustBeInteger,mustBeGreaterThan(options.LookupPoints,100)} = 801
    options.IntegrationPoints (1,1) double {mustBeInteger,mustBeGreaterThan(options.IntegrationPoints,200)} = 1200
    options.FoundationModel (1,1) string = "none"
    options.GlobalScale (1,1) double {mustBePositive} = 1
end

foundationModel = lower(strtrim(options.FoundationModel));
if abs(options.GlobalScale-1) > 10*eps
    error('pg_build_pengyue_pr_tvms_model:GlobalScaleForbidden', ...
        'GlobalScale is disabled for the physical reproduction route; use 1.');
end
if foundationModel == "constant"
    error('pg_build_pengyue_pr_tvms_model:ConstantFoundationForbidden', ...
        ['The temporary constant foundation compliance is disabled. ' ...
         'Use FoundationModel="none" or provide explicit ring-body geometry.']);
end
if ~ismember(foundationModel,["none","chaari"])
    error('pg_build_pengyue_pr_tvms_model:FoundationModel', ...
        'FoundationModel must be "none" or "chaari".');
end
if foundationModel == "chaari"
    error('pg_build_pengyue_pr_tvms_model:RingFoundationGeometryMissing', ...
        ['A ring-specific Chaari/Sainsot foundation model cannot be built ' ...
         'until rim/body and tooth-root geometry are supplied.']);
end

zs = p.gear.zs; %#ok<NASGU>
zp = p.gear.zp;
zr = p.gear.zr;
m = p.gear.module;
alpha0 = p.gear.alpha_nominal;
E = p.material.gear.E;
nu = p.material.gear.nu;
G = p.material.gear.G;
Bp = p.gear.faceWidthSunPlanet;
Br = p.gear.faceWidthRing;

% Reported/derived standard geometry.
rp = p.gear.pitchRadiusPlanet;
rr = p.gear.pitchRadiusRing;
rbP = p.gear.baseRadiusPlanet;
rbR = p.gear.baseRadiusRing;
raP = p.gear.addendumRadiusPlanet;
raR = p.gear.addendumRadiusRing;  % internal-gear addendum is the inner tip
rfP = p.gear.rootRadiusPlanet;
rfR = p.gear.rootRadiusRing;      % internal-gear dedendum is the outer root
aRP = rr-rp;
basePitch = pi*m*cos(alpha0);

% Internal-mesh path of contact (same sign convention as the thesis).
pathApproach = sqrt(max(raP^2-rbP^2,0)) - sqrt(max(raR^2-rbR^2,0));
pathRecess = aRP*sin(alpha0);
pathTotal = pathApproach + pathRecess;
epsilon = pathTotal/basePitch;
if ~(isfinite(epsilon) && epsilon > 1 && epsilon < 2)
    error('pg_build_pengyue_pr_tvms_model:ContactRatio', ...
        'Internal-mesh contact ratio %.9g is not in (1,2).',epsilon);
end

geom = resolveInternalGeometry(p,alpha0,zr,rbR,raR,rfR);

eta = linspace(0,1,options.LookupPoints).';
kHealthy = zeros(size(eta));
kFault = zeros(size(eta));
% [Hertz, planet bending, planet axial, planet shear, planet foundation,
%  ring bending, ring axial, ring shear, ring foundation, total]
compHealthy = zeros(numel(eta),10);
compFault = zeros(numel(eta),10);
ringDetail = repmat(emptyRingDetail(),numel(eta),1);

for i = 1:numel(eta)
    % Map the normalized contact coordinate onto the two involute contact
    % radii.  For an internal gear the addendum circle is the inner circle,
    % so the ring tangent length increases from r_a,r to r_f,r.  Using the
    % tangent lengths explicitly avoids accidentally clamping the whole
    % ring path to r_a,r when the line-of-action origin is changed.
    lPa = sqrt(max(raP^2-rbP^2,0));
    lPf = sqrt(max(rfP^2-rbP^2,0));
    lRa = sqrt(max(raR^2-rbR^2,0));
    lRf = sqrt(max(rfR^2-rbR^2,0));
    lP = lPa + eta(i)*(lPf-lPa);
    lR = lRa + eta(i)*(lRf-lRa);
    rCP = sqrt(rbP^2+lP^2);
    rCR = sqrt(rbR^2+lR^2);

    plaH = externalPlanetCompliance(rCP,p,options.IntegrationPoints);
    % A planet-root crack, when requested, modifies only planet bending and
    % shear exactly as in the SP model.  A sun crack leaves PR healthy.
    plaC = externalPlanetCompliance(rCP,p,options.IntegrationPoints, ...
        isPlanetFault(fault),fault);
    [ringC,rd] = internalRingCompliance(rCR,p,geom,options.IntegrationPoints);
    ringDetail(i) = rd;

    kHertz = pi*E*Br/(4*(1-nu^2));
    cH = 1/kHertz;
    cHealthy = cH + plaH.bending + plaH.axial + plaH.shear + 0 + ...
        ringC.bending + ringC.axial + ringC.shear + 0;
    cFault = cH + plaC.bending + plaH.axial + plaC.shear + 0 + ...
        ringC.bending + ringC.axial + ringC.shear + 0;

    kHealthy(i) = options.GlobalScale/cHealthy;
    kFault(i) = options.GlobalScale/cFault;
    compHealthy(i,:) = [cH plaH.bending plaH.axial plaH.shear 0 ...
        ringC.bending ringC.axial ringC.shear 0 cHealthy];
    compFault(i,:) = [cH plaC.bending plaH.axial plaC.shear 0 ...
        ringC.bending ringC.axial ringC.shear 0 cFault];
end

planetPhase = p.gear.planetPhase0(:).';
meshPhaseCycles = mod(zr*planetPhase/(2*pi),1);
toothPhaseCycles = zr*planetPhase/(2*pi);

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
model.foundationModel = foundationModel;
model.ringFoundationStatus = 'MISSING: ring rim/body and tooth-root geometry; compliance set to 0';
model.internalGeometry = geom;
model.ringDetail = ringDetail;
model.complianceColumnNames = {'hertz','planetBending','planetAxial','planetShear', ...
    'planetFoundation','ringBending','ringAxial','ringShear','ringFoundation','total'};
model.complianceHealthy = compHealthy;
model.complianceFault = compFault;
model.method = ['Potential-energy internal-mesh TVMS. Planet tooth uses ' ...
    'the SP numerical section integration; ring tooth uses Peng Yue ' ...
    'Eqs. (2-25)-(2-27). Foundation term is disabled because ring-body ' ...
    'geometry is missing.'];
end

function tf = isPlanetFault(fault)
tf = isfield(fault,'type') && strcmpi(fault.type,'planet') && ...
    isfield(fault,'crackDepth') && fault.crackDepth > 0;
end

function c = externalPlanetCompliance(rF,p,nInt,applyCrack,fault)
% Numerical equivalent of the external tooth potential-energy terms.
if nargin < 4; applyCrack = false; end
if nargin < 5; fault = struct(); end
rb = p.gear.baseRadiusPlanet;
rf = p.gear.rootRadiusPlanet;
z = p.gear.zp;
alpha0 = p.gear.alpha_nominal;
E = p.material.gear.E; G = p.material.gear.G; B = p.gear.faceWidthSunPlanet;
phiB = pi/(2*z) + invInvolute(alpha0);
r = linspace(rf,rF,nInt).';
alphaR = zeros(size(r));
mask = r > rb;
alphaR(mask) = acos(rb./r(mask));
thetaHalf = phiB - invInvolute(alphaR);
thetaHalf(~mask) = phiB;
xAxis = r.*cos(thetaHalf);
xRoot = rf*cos(phiB);
x = xAxis-xRoot;
hx = r.*sin(thetaHalf);
sHealthy = max(2*hx,1e-9);
sSection = sHealthy;
if applyCrack
    q = fault.crackDepth;
    gamma = fault.crackAngle;
    hc = rf*sin(phiB);
    hq = hc-q*sin(gamma);
    affected = hx > hq;
    sCrack = hx+hc-q*sin(gamma);
    sSection(affected) = max(sCrack(affected),0.005*sHealthy(affected));
    sSection = min(sSection,sHealthy);
end
A = B*sHealthy;
AS = B*sSection;
IS = B*sSection.^3/12;
alphaF = acos(min(max(rb/rF,0),1));
xF = x(end); hF = hx(end);
Fa = sin(alphaF); Fb = cos(alphaF);
M = Fb.*(xF-x)-Fa.*hF;
[xu,ia] = unique(x,'stable');
M=M(ia); A=A(ia); AS=AS(ia); IS=IS(ia);
c.bending=trapz(xu,M.^2./(E*IS));
c.axial=trapz(xu,Fa^2./(E*A));
c.shear=trapz(xu,1.2*Fb^2./(G*AS));
end

function [c,d] = internalRingCompliance(rContact,p,g,nInt)
% Peng Yue Eqs. (2-25)-(2-27), evaluated for the current ring contact.
E=p.material.gear.E; G=p.material.gear.G; B=p.gear.faceWidthRing;
rb=p.gear.baseRadiusRing;
alphaF=acos(min(max(rb/rContact,0),1));
alphaLower=g.alphaRoot;
alphaUpper=alphaF;
% The internal coordinate is oppositely oriented to the external one; use
% an ordered integration interval and retain the positive energy integral.
lo=min(alphaLower,alphaUpper); hi=max(alphaLower,alphaUpper);
if hi-lo < 10*eps
    c=struct('bending',0,'axial',0,'shear',0);
    denMin=NaN;
else
    ax=linspace(lo,hi,nInt).';
    D1=sin(g.theta+ax); D2=cos(g.theta+ax);
    D3=sin(g.theta+alphaF); D4=cos(g.theta+alphaF);
    D=D4.*(D2+(g.phiB+ax).*D1-D4-(g.phiB+alphaF).*D3) ...
      -D3.*(D3-(g.phiB+alphaF).*D4);
    den=2*D1-2*(ax+g.phiB).*D2;
    denSafe=den;
    small=abs(denSafe)<1e-12;
    denSafe(small)=1e-12;
    ib= D.^2.*12.*(-ax-g.phiB).*D2./denSafe.^3;
    ia= D3.^2.*(-ax-g.phiB).*D2./denSafe;
    is= 1.2*D4.^2.*(-alphaF-g.phiB).*D2./denSafe;
    c=struct('bending',abs(trapz(ax,ib)/(E*B)), ...
        'axial',abs(trapz(ax,ia)/(E*B)), ...
        'shear',abs(trapz(ax,is)/(G*B)));
    denMin=min(abs(den));
    if any(~isfinite([c.bending c.axial c.shear]))
        error('pg_build_pengyue_pr_tvms_model:InternalCompliance', ...
            'Non-finite ring compliance at rContact=%.9g m.',rContact);
    end
end
d=struct('rContact',rContact,'baseRadius',rb,'alphaF',alphaF, ...
    'alphaRoot',g.alphaRoot,'theta',g.theta,'phiB',g.phiB, ...
    'denominatorMin',denMin,'formula','Peng Yue Eqs. (2-25)-(2-27)', ...
    'foundationStatus','MISSING: ring foundation geometry');
end

function g=resolveInternalGeometry(p,alpha0,zr,rbR,raR,rfR)
g=struct();
g.geometryStatus='inferred standard internal involute (thesis does not tabulate angles)';
g.theta=alpha0;
g.phiB=pi/(2*zr)-invInvolute(alpha0);
g.alphaRoot=acos(min(max(rbR/rfR,0),1));
if isfield(p.gear,'internalMeshGeometry')
    q=p.gear.internalMeshGeometry;
    if isfield(q,'theta') && isfinite(q.theta); g.theta=q.theta; end
    if isfield(q,'phiB') && isfinite(q.phiB); g.phiB=q.phiB; end
    if isfield(q,'alphaRoot') && isfinite(q.alphaRoot); g.alphaRoot=q.alphaRoot; end
    g.geometryStatus='explicit p.gear.internalMeshGeometry override';
end
g.rb=rbR; g.ra=raR; g.rf=rfR; g.alphaAddendum=acos(min(max(rbR/raR,0),1));
end

function y=invInvolute(a)
y=tan(a)-a;
end

function d=emptyRingDetail()
d=struct('rContact',NaN,'baseRadius',NaN,'alphaF',NaN,'alphaRoot',NaN, ...
    'theta',NaN,'phiB',NaN,'denominatorMin',NaN,'formula','', ...
    'foundationStatus','');
end
