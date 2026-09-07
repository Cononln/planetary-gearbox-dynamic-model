function truth = pg_fault_event_truth(time,stiffnessLoss,p,fault)
%PG_FAULT_EVENT_TRUTH Extract source-event truth from a simulated defect.
% Each complete nonzero stiffness-loss interval is one known fault-contact
% event.  Its occurrence time is the maximum stiffness-loss instant, rather
% than an arbitrary cycle origin.  The returned source time is independent
% of any sensor response and can therefore be used to score event-arrival
% recovery objectively.

time = time(:);
if isempty(time) || any(~isfinite(time)) || any(diff(time)<=0)
    error('pg_fault_event_truth:Time', ...
        'time must be a finite strictly increasing vector.');
end
if isvector(stiffnessLoss)
    stiffnessLoss = stiffnessLoss(:);
end
if size(stiffnessLoss,1)~=numel(time) || any(~isfinite(stiffnessLoss(:)))
    error('pg_fault_event_truth:Loss', ...
        'stiffnessLoss must be finite with one row per time sample.');
end
if any(stiffnessLoss(:)<-eps)
    error('pg_fault_event_truth:LossSign', ...
        'stiffnessLoss must be nonnegative.');
end

combinedLoss = sum(max(stiffnessLoss,0),2);
threshold = max(combinedLoss)*1e-10;
active = combinedLoss>max(threshold,eps(max(combinedLoss)));
edge = diff([false;active;false]);
startIndex = find(edge==1);
endIndex = find(edge==-1)-1;
complete = startIndex>1 & endIndex<numel(time);
startIndexComplete = startIndex(complete);
endIndexComplete = endIndex(complete);
nEvent = numel(startIndexComplete);

eventIndex = zeros(nEvent,1);
planetIndex = zeros(nEvent,1);
peakLoss = zeros(nEvent,1);
for i = 1:nEvent
    interval = startIndexComplete(i):endIndexComplete(i);
    [peakLoss(i),localIndex] = max(combinedLoss(interval));
    eventIndex(i) = interval(localIndex);
    [~,planetIndex(i)] = max(stiffnessLoss(eventIndex(i),:));
end
motion = pg_motion_state(time(eventIndex),p);

truth.eventTime = time(eventIndex);
truth.sampleIndex = eventIndex;
truth.carrierAngleRad = motion.phi_c(:);
truth.carrierAngleDeg = mod(rad2deg(motion.phi_c(:)),360);
truth.meshPhaseRad = motion.meshPhaseRad(:);
truth.planetIndex = planetIndex;
truth.peakStiffnessLossNPerM = peakLoss;
truth.eventStartTime = time(startIndexComplete);
truth.eventEndTime = time(endIndexComplete);
truth.completeEventCount = nEvent;
truth.clippedEventCount = nnz(~complete);
truth.activeMask = active;
truth.fault = fault;
truth.definition = [ ...
    'Source event time is the maximum simulated defect stiffness-loss ' ...
    'instant within each complete fault-contact interval.'];
end
