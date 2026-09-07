function state = pg_motion_state(t,p)
%PG_MOTION_STATE Kinematically consistent prescribed carrier motion.
% The optional speed profile is an externally prescribed operating
% trajectory.  It drives every moving-coordinate quantity in the model:
% carrier position, mesh counter, mesh frequency, TVMS, and transmission
% error.  It is intentionally separate from the small torque ripple used
% in the force model.

if ~isnumeric(t) || any(~isfinite(t(:)))
    error('pg_motion_state:Time','t must contain finite numeric values.');
end

tShape = size(t);
t = double(t(:));
phi0 = p.model.phi_c0;
omegaCarrierNominal = p.kin.omega_c;

enabled = false;
fraction = 0;
frequencyHz = 0;
phaseRad = 0;
if isfield(p,'operating') && isfield(p.operating,'speedProfile')
    profile = p.operating.speedProfile;
    if isfield(profile,'enabled'); enabled = logical(profile.enabled); end
    if isfield(profile,'fraction'); fraction = profile.fraction; end
    if isfield(profile,'frequencyHz'); frequencyHz = profile.frequencyHz; end
    if isfield(profile,'phaseRad'); phaseRad = profile.phaseRad; end
end

if ~isscalar(fraction) || ~isfinite(fraction) || abs(fraction)>=1
    error('pg_motion_state:Fraction', ...
        'speedProfile.fraction must be finite and have magnitude below one.');
end
if ~isscalar(frequencyHz) || ~isfinite(frequencyHz) || frequencyHz<0
    error('pg_motion_state:Frequency', ...
        'speedProfile.frequencyHz must be a finite nonnegative scalar.');
end
if ~isscalar(phaseRad) || ~isfinite(phaseRad)
    error('pg_motion_state:Phase', ...
        'speedProfile.phaseRad must be a finite scalar.');
end

if enabled && fraction~=0
    if frequencyHz>0
        modulationPhase = 2*pi*frequencyHz*t+phaseRad;
        omegaCarrier = omegaCarrierNominal*(1+fraction*cos(modulationPhase));
        phiCarrier = phi0+omegaCarrierNominal*(t+ ...
            fraction*(sin(modulationPhase)-sin(phaseRad))/(2*pi*frequencyHz));
    else
        omegaCarrier = omegaCarrierNominal*(1+fraction*cos(phaseRad));
        phiCarrier = phi0+omegaCarrier*t;
    end
else
    omegaCarrier = omegaCarrierNominal+zeros(size(t));
    phiCarrier = phi0+omegaCarrierNominal*t;
end

% For a fixed ring, omega_mesh = Zr*omega_c and
% omega_s = (Zs+Zr)/Zs*omega_c.  The mesh phase is measured from t = 0 so
% the existing tooth-indexing phases remain valid.
omegaSun = (p.gear.zs+p.gear.zr)/p.gear.zs*omegaCarrier;
omegaMesh = p.gear.zr*omegaCarrier;
meshPhaseRad = p.gear.zr*(phiCarrier-phi0);

state.phi_c = reshape(phiCarrier,tShape);
state.omega_c = reshape(omegaCarrier,tShape);
state.omega_s = reshape(omegaSun,tShape);
state.omega_mesh = reshape(omegaMesh,tShape);
state.f_c = state.omega_c/(2*pi);
state.f_s = state.omega_s/(2*pi);
state.f_mesh = state.omega_mesh/(2*pi);
state.meshPhaseRad = reshape(meshPhaseRad,tShape);
state.meshCounter = state.meshPhaseRad/(2*pi);
state.speedProfileEnabled = enabled;
state.speedProfileFraction = fraction;
state.speedProfileFrequencyHz = frequencyHz;
end
