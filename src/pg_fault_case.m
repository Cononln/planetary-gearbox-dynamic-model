function fault = pg_fault_case(type, severity, p)
%PG_FAULT_CASE Create a frozen partial-width broken-tooth definition.
% type is none, sun, or planet. Severity is the missing face-width fraction.
% Use pg_root_crack_case for the energy-based sun-root crack model.

if nargin < 3
    p = pg_parameters;
end
type = lower(string(type));
validType = any(type == ["none","sun","planet"]);
if ~validType
    error('pg_fault_case:Type','Fault type must be none, sun, or planet.');
end
if ~isscalar(severity) || severity < 0 || severity > 1
    error('pg_fault_case:Severity','Severity must lie in [0,1].');
end
if type == "none" && severity ~= 0
    error('pg_fault_case:HealthySeverity','Healthy case must have severity 0.');
end

fault.type = char(type);
if type == "none"
    fault.damageModel = 'none';
else
    fault.damageModel = 'partial_width_break';
end
fault.severity = severity;
fault.tooth = p.tvms.faultTooth;
fault.targetPlanet = p.tvms.targetPlanet;
fault.fullFaceWidth = p.gear.faceWidthSunPlanet;
fault.missingFaceWidth = severity*fault.fullFaceWidth;
fault.remainingFaceWidth = (1-severity)*fault.fullFaceWidth;
fault.radialRemovalDepth = p.tvms.radialRemovalDepth;
fault.depthDefinition = 'removed to root circle';
fault.speedRippleFraction = 0;
fault.speedRipplePhase = 0;
end
