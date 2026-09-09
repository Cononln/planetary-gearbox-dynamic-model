function audit = tvms_geometry_audit(c)
%TVMS_GEOMETRY_AUDIT Check planetary closure and required TVMS geometry.
tol = 1e-9;
audit.standardRingTeeth = c.zs + 2*c.zp;
audit.standardExternalCenter = c.module*(c.zs+c.zp)/2;
audit.internalCenter = c.module*(c.zr-c.zp)/2;
audit.externalCenterResidual = c.centerSP-audit.standardExternalCenter;
audit.internalCenterResidual = c.centerPR-audit.internalCenter;
% `Zr=Zs+2*Zp` is only the zero-centre-modification special case.  For a
% shared planet centre, the two centre modifications close the tooth count:
% Zr = Zs + 2*Zp + 2*(ySP-yPR), y=(a-a0)/m.
audit.centerModificationSP = audit.externalCenterResidual/c.module;
audit.centerModificationPR = audit.internalCenterResidual/c.module;
audit.centerModifiedRingTeeth = audit.standardRingTeeth + ...
    2*(audit.centerModificationSP-audit.centerModificationPR);
audit.toothClosureResidual = c.zr-audit.centerModifiedRingTeeth;
audit.standardZeroModification = abs(audit.centerModificationSP)<tol && ...
    abs(audit.centerModificationPR)<tol;
audit.centresCloseToothCount = abs(audit.toothClosureResidual)<1e-6;

% A stated working centre constrains a profile-shift combination under the
% standard involute tooth-thickness convention.  A mismatch can also mean
% unknown backlash, tooth-thickness correction, or non-standard tooth form.
invAlpha = @(alpha) tan(alpha)-alpha;
audit.alphaWorkingSP = acos(audit.standardExternalCenter * ...
    cos(c.alphaNominal)/c.centerSP);
audit.impliedProfileShiftSumSP = (c.zs+c.zp)/(2*tan(c.alphaNominal)) * ...
    (invAlpha(audit.alphaWorkingSP)-invAlpha(c.alphaNominal));
if isfinite(c.xSun) && isfinite(c.xPlanet)
    audit.suppliedProfileShiftSumSP = c.xSun+c.xPlanet;
    audit.profileShiftResidualSP = audit.suppliedProfileShiftSumSP - ...
        audit.impliedProfileShiftSumSP;
    audit.profileShiftSPConsistent = abs(audit.profileShiftResidualSP)<1e-3;
else
    audit.suppliedProfileShiftSumSP = NaN;
    audit.profileShiftResidualSP = NaN;
    audit.profileShiftSPConsistent = true;
end
audit.alphaWorkingPR = acos(audit.internalCenter * ...
    cos(c.alphaNominal)/c.centerPR);
audit.impliedProfileShiftDifferencePR = (c.zr-c.zp)/(2*tan(c.alphaNominal)) * ...
    (invAlpha(audit.alphaWorkingPR)-invAlpha(c.alphaNominal));
if isfinite(c.xRing) && isfinite(c.xPlanet)
    audit.suppliedProfileShiftDifferencePR = c.xRing-c.xPlanet;
    audit.profileShiftResidualPR = audit.suppliedProfileShiftDifferencePR - ...
        audit.impliedProfileShiftDifferencePR;
    audit.profileShiftPRConsistent = abs(audit.profileShiftResidualPR)<1e-3;
else
    audit.suppliedProfileShiftDifferencePR = NaN;
    audit.profileShiftResidualPR = NaN;
    audit.profileShiftPRConsistent = true;
end
missing = {'ring.tipRadius','ring.rootRadius','haStar','clearanceStar', ...
    'xSun','xPlanet','xRing','foundationCalibration'};
values = [c.ring.tipRadius,c.ring.rootRadius,c.haStar,c.clearanceStar, ...
    c.xSun,c.xPlanet,c.xRing,c.foundationCalibration];
audit.missing = missing(~isfinite(values));
% The present project is kinematically closed by its specified carrier
% radius.  The missing item is the *allocation* of the implied external
% profile shift, plus the internal-tooth dimensions needed by the energy
% integrals; it is not a reason to alter the inherited tooth counts.
audit.modelGeometryExplained = audit.centresCloseToothCount;
audit.profileGeometryConsistent = audit.profileShiftSPConsistent && ...
    audit.profileShiftPRConsistent;
if isfield(c,'useEquivalentWorkingCentreThicknessCorrection') && ...
        c.useEquivalentWorkingCentreThicknessCorrection
    audit.profileGeometryConsistent = audit.profileShiftPRConsistent;
    audit.equivalentThicknessCorrectionSP = audit.impliedProfileShiftSumSP;
else
    audit.equivalentThicknessCorrectionSP = NaN;
end
audit.ready = isempty(audit.missing) && audit.modelGeometryExplained && ...
    audit.profileGeometryConsistent;
if ~audit.centresCloseToothCount
    audit.warning = sprintf(['WARNING: tooth counts and stated centres do not close: ' ...
        'Zr=%g, centre-modified closure=%g.'], ...
        c.zr,audit.centerModifiedRingTeeth);
elseif ~audit.profileGeometryConsistent
    audit.warning = sprintf(['Provided profile shifts do not account for the stated ' ...
        'working centres under standard tooth-thickness convention: SP supplied ' ...
        'xsum=%g versus centre-implied %g.  Obtain tooth thickness/backlash or ' ...
        'tooth-form correction data before TVMS.'], ...
        audit.suppliedProfileShiftSumSP,audit.impliedProfileShiftSumSP);
elseif ~audit.standardZeroModification
    audit.warning = sprintf(['Modified-centre planetary geometry: zero-modification reference ' ...
        'Zr=%g, centre-modified closure=%g.  Profile shifts, tooth thickness, ' ...
        'backlash and tooth-form data must be supplied before TVMS is computed.'], ...
        audit.standardRingTeeth,audit.centerModifiedRingTeeth);
else
    audit.warning = '';
end
end
