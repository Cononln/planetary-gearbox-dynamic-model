function aligned90 = pg_apply_phase_truth(spectrum90,truth)
%PG_APPLY_PHASE_TRUTH Phase-align an angle-frequency spectrum to sensor 0.
% spectrum90 must have the same frequency-by-carrier-angle layout as truth.

if ~isequal(size(spectrum90),size(truth.phaseDifferenceWrapped))
    error('pg_apply_phase_truth:Size', ...
        'Input spectrum must match the truth frequency-by-angle grid.');
end
aligned90 = spectrum90.*exp(-1i*truth.phaseDifferenceWrapped);
end
