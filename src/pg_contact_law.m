function [force,tangentK,tangentC,isClosed] = pg_contact_law( ...
    delta,deltaDot,kMesh,cMesh,model,clearance,dampingOnlyInContact)
%PG_CONTACT_LAW Mesh force for linear, unilateral, or backlash contact.
% Positive delta denotes compression. tangentK and tangentC are positive
% derivatives with respect to compression and compression rate.

if nargin < 7
    dampingOnlyInContact = true;
end
if clearance < 0 || ~isfinite(clearance)
    error('pg_contact_law:Clearance','clearance must be finite and >= 0.');
end
model = lower(string(model));
switch model
    case "linear"
        elasticDeflection = delta;
        tangentK = ones(size(delta));
        isClosed = true(size(delta));
    case "unilateral"
        isClosed = delta > clearance;
        elasticDeflection = max(delta-clearance,0);
        tangentK = double(isClosed);
    case "backlash"
        positiveFlank = delta > clearance;
        negativeFlank = delta < -clearance;
        isClosed = positiveFlank | negativeFlank;
        elasticDeflection = zeros(size(delta));
        elasticDeflection(positiveFlank) = delta(positiveFlank)-clearance;
        elasticDeflection(negativeFlank) = delta(negativeFlank)+clearance;
        tangentK = double(isClosed);
    otherwise
        error('pg_contact_law:Model','Unsupported contact model: %s',model);
end

if dampingOnlyInContact
    tangentC = double(isClosed);
else
    tangentC = ones(size(delta));
end
force = kMesh.*elasticDeflection+cMesh.*tangentC.*deltaDot;
tangentK = kMesh.*tangentK;
tangentC = cMesh.*tangentC;
end
