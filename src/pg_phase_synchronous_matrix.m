function out = pg_phase_synchronous_matrix(signal,phase,samplesPerCycle)
%PG_PHASE_SYNCHRONOUS_MATRIX Resample a signal into phase-locked cycles.
% phase is any recovered monotonically increasing phase, e.g. mesh phase or
% carrier phase.  This is the interface used by TSA/TPSVD; it needs no
% encoder once phase has been estimated.

if isvector(signal)
    signal = signal(:);
end
phase = phase(:);
phase = phase-phase(1);
firstCycle = ceil(phase(1)/(2*pi));
lastCycle = floor(phase(end)/(2*pi))-1;
cycleIndex = (firstCycle:lastCycle).';
withinCycle = (0:samplesPerCycle-1)/samplesPerCycle;
queryPhase = 2*pi*(cycleIndex+withinCycle);
nCycle = size(queryPhase,1);
nSensor = size(signal,2);
matrix = zeros(nCycle,samplesPerCycle,nSensor);
for is = 1:nSensor
    values = interp1(phase,signal(:,is),queryPhase(:),'linear');
    matrix(:,:,is) = reshape(values,size(queryPhase));
end

out.matrix = matrix;
out.cycleIndex = cycleIndex;
out.withinCyclePhase = 2*pi*withinCycle;
out.queryPhase = queryPhase;
out.samplesPerCycle = samplesPerCycle;
end
