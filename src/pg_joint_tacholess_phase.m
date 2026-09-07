function result = pg_joint_tacholess_phase(component,fs,harmonicOrder, ...
    nominalMeshFrequencyHz,ringToothCount,options)
%PG_JOINT_TACHOLESS_PHASE Joint common/path phase separation.
%
% phase(z_s,h) = h*phi_m + psi_s,h(phi_m/Zr) + error
%
% phi_m is common to all synchronized sensors and harmonics.  Each path
% phasor is fitted in carrier-angle spatial orders that are integer
% multiples of options.spatialMultiple (three for three equally spaced
% planets).  With three or more sensors, every sensor pair is fitted and
% projected onto a zero-mean graph of sensor path phases.  The graph
% projection enforces pairwise cycle consistency and provides a closure
% residual for reliability assessment.  The fitted phasor is divided out
% without changing amplitude.  Absolute constant phase is deliberately
% left unidentified.

if nargin < 6
    options = struct;
end
options = fillDefaults(options);
harmonicOrder = harmonicOrder(:).';
[nSample,nSensor,nHarmonic] = size(component);
if nHarmonic ~= numel(harmonicOrder)
    error('pg_joint_tacholess_phase:HarmonicSize', ...
        'The third component dimension must match harmonicOrder.');
end
if nSensor < 2
    error('pg_joint_tacholess_phase:SensorCount', ...
        'At least two synchronized vibration channels are required.');
end

t = (0:nSample-1)'/fs;
nominalPhase = 2*pi*nominalMeshFrequencyHz*t;
amplitudeWeight = reliabilityWeights(component);
commonWeight = amplitudeWeight;
pathWeight = amplitudeWeight;
if ~isempty(options.sensorHarmonicWeight)
    sensorHarmonicWeight = options.sensorHarmonicWeight;
    if ~isequal(size(sensorHarmonicWeight),[nSensor,nHarmonic])
        error('pg_joint_tacholess_phase:ReliabilitySize', ...
            'sensorHarmonicWeight must be nSensor-by-nHarmonic.');
    end
    commonWeight = commonWeight.*reshape( ...
        sensorHarmonicWeight,[1,nSensor,nHarmonic]);
    pathWeight = pathWeight.*reshape( ...
        sensorHarmonicWeight,[1,nSensor,nHarmonic]);
end
if ~isempty(options.commonSensorHarmonicWeight)
    commonSensorHarmonicWeight = options.commonSensorHarmonicWeight;
    if ~isequal(size(commonSensorHarmonicWeight),[nSensor,nHarmonic])
        error('pg_joint_tacholess_phase:CommonReliabilitySize', ...
            'commonSensorHarmonicWeight must be nSensor-by-nHarmonic.');
    end
    commonWeight = commonWeight.*reshape( ...
        commonSensorHarmonicWeight,[1,nSensor,nHarmonic]);
end
if ~isempty(options.pathSensorHarmonicWeight)
    pathSensorHarmonicWeight = options.pathSensorHarmonicWeight;
    if ~isequal(size(pathSensorHarmonicWeight),[nSensor,nHarmonic])
        error('pg_joint_tacholess_phase:PathReliabilitySize', ...
            'pathSensorHarmonicWeight must be nSensor-by-nHarmonic.');
    end
    pathWeight = pathWeight.*reshape( ...
        pathSensorHarmonicWeight,[1,nSensor,nHarmonic]);
end

% Estimate one-sample phase increments around the nominal trajectory.  This
% avoids independently unwrapping six absolute phases and then inheriting a
% permanent cycle slip from one low-amplitude path anti-resonance.
phaseEstimate = incrementalCommonPhase(component,harmonicOrder, ...
    nominalPhase,commonWeight,fs);
phaseEstimate = regularizeCommonPhase(phaseEstimate,nominalPhase, ...
    phaseEstimate/ringToothCount,fs,options);
phaseEstimate = phaseEstimate+2*pi*options.meanFrequencyCorrectionHz*t;
phaseAnchor = phaseEstimate;

history = zeros(nSample,options.nIteration+1);
history(:,1) = phaseEstimate;
if isempty(options.spatialOrderOverride)
    spatialOrder = options.spatialMultiple* ...
        (-options.pathOrder:options.pathOrder);
else
    spatialOrder = options.spatialOrderOverride(:).';
end
pathPenalty = ones(size(spatialOrder));
isOffGrid = mod(abs(spatialOrder),options.spatialMultiple)~=0;
pathPenalty(isOffGrid) = options.offGridPenalty;

for iteration = 1:options.nIteration
    carrierPhase = phaseEstimate/ringToothCount;
    basis = exp(1i*carrierPhase*spatialOrder);
    [~,corrected,~,~] = fitDifferentialPaths( ...
        component,basis,pathWeight,options.ridge,pathPenalty, ...
        options.minEdgeFitCoherence);

    update = incrementalCommonPhase(corrected,harmonicOrder, ...
        phaseEstimate,commonWeight,fs);
    update = regularizeCommonPhase(update,nominalPhase,carrierPhase, ...
        fs,options);
    tentative = (1-options.relaxation)*phaseEstimate+ ...
        options.relaxation*update;
    % The raw multiharmonic common mode fixes the remaining common-path
    % gauge.  Differential-path fitting is allowed to refine, not replace,
    % that robust speed estimate.
    phaseEstimate = options.commonAnchorWeight*phaseAnchor+ ...
        (1-options.commonAnchorWeight)*tentative;
    history(:,iteration+1) = phaseEstimate;
end

% Refit once at the returned common phase so paths and phase are consistent.
carrierPhase = phaseEstimate/ringToothCount;
basis = exp(1i*carrierPhase*spatialOrder);
[pathPhasor,corrected,coefficient,pathGraph] = fitDifferentialPaths( ...
    component,basis,pathWeight,options.ridge,pathPenalty, ...
    options.minEdgeFitCoherence);

instantaneousMeshFrequencyHz = [diff(phaseEstimate); ...
    phaseEstimate(end)-phaseEstimate(end-1)]*fs/(2*pi);
instantaneousMeshFrequencyHz = movingAverage( ...
    instantaneousMeshFrequencyHz,max(3,round(options.frequencySmoothSec*fs)));

result.meshPhase = phaseEstimate;
result.carrierPhase = carrierPhase;
result.instantaneousMeshFrequencyHz = instantaneousMeshFrequencyHz;
result.pathPhasor = pathPhasor;
result.correctedComponent = corrected;
result.correctedSignal = squeeze(real(sum(corrected,3)));
result.amplitudeWeight = amplitudeWeight;
result.commonWeight = commonWeight;
result.pathWeight = pathWeight;
result.spatialOrder = spatialOrder;
result.pathPenalty = pathPenalty;
result.pathCoefficient = coefficient;
result.pathGraph = pathGraph;
result.iterationHistory = history;
result.options = options;
result.identifiability = [ ...
    'Phase is identifiable up to a constant. Differential moving-path ', ...
    'phase is identifiable from synchronized channels; a path phase that ', ...
    'is exactly common to every channel remains in the common-mode gauge.'];
end

function options = fillDefaults(options)
defaults.nIteration = 3;
defaults.spatialMultiple = 3;
defaults.pathOrder = 4;
defaults.ridge = 1e-3;
defaults.relaxation = 0.35;
defaults.commonAnchorWeight = 0.80;
defaults.phaseSmoothSec = 0.012;
defaults.frequencySmoothSec = 0.040;
defaults.removeCommonSpatialOrder = false;
defaults.spatialOrderOverride = [];
defaults.offGridPenalty = 1;
defaults.meanFrequencyCorrectionHz = 0;
defaults.sensorHarmonicWeight = [];
defaults.commonSensorHarmonicWeight = [];
defaults.pathSensorHarmonicWeight = [];
defaults.minEdgeFitCoherence = 0;
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(options,names{i}) || isempty(options.(names{i}))
        options.(names{i}) = defaults.(names{i});
    end
end
end

function [pathPhasor,corrected,coefficient,graph] = fitDifferentialPaths( ...
    component,basis,amplitudeWeight,ridge,pathPenalty,minEdgeFitCoherence)
[nSample,nSensor,nHarmonic] = size(component);
pathPhasor = complex(ones(size(component)));
corrected = component;
[pairIndex,incidence] = completePairGraph(nSensor);
nPair = size(pairIndex,1);
coefficient = cell(nPair,nHarmonic);
graph.pairIndex = pairIndex;
graph.incidence = incidence;
graph.edgeFitRmseDeg = zeros(nPair,nHarmonic);
graph.edgeFitCoherence = zeros(nPair,nHarmonic);
graph.edgeWeight = zeros(nPair,nHarmonic);
graph.activeEdgeMask = false(nPair,nHarmonic);
graph.activeSensorMask = false(nSensor,nHarmonic);
graph.edgeProjectionRmseDeg = zeros(nHarmonic,1);
graph.cycleClosureRmseDeg = nan(nHarmonic,1);
for ih = 1:nHarmonic
    edgePhase = zeros(nSample,nPair);
    edgeWeight = zeros(nPair,1);
    for ip = 1:nPair
        sensorI = pairIndex(ip,1);
        sensorJ = pairIndex(ip,2);
        crossPhasor = component(:,sensorJ,ih).* ...
            conj(component(:,sensorI,ih));
        crossPhasor = crossPhasor./max(abs(crossPhasor),eps);
        weight = sqrt(amplitudeWeight(:,sensorJ,ih).* ...
            amplitudeWeight(:,sensorI,ih));
        if sum(weight)<=eps
            coefficient{ip,ih} = complex(zeros(size(basis,2),1));
            coefficient{ip,ih}(ceil(size(basis,2)/2)) = 1;
            edgePhase(:,ip) = 0;
            graph.edgeFitRmseDeg(ip,ih) = nan;
            graph.edgeFitCoherence(ip,ih) = 0;
            edgeWeight(ip) = eps;
            continue;
        end
        coefficient{ip,ih} = weightedComplexFit( ...
            basis,crossPhasor,weight,ridge,pathPenalty);
        fitted = basis*coefficient{ip,ih};
        fitted = fitted./max(abs(fitted),eps);
        edgePhase(:,ip) = unwrap(angle(fitted));
        fitError = angle(crossPhasor.*conj(fitted));
        fitRmse = sqrt(sum(weight.*fitError.^2)/max(sum(weight),eps));
        graph.edgeFitRmseDeg(ip,ih) = rad2deg(fitRmse);
        graph.edgeFitCoherence(ip,ih) = abs(sum( ...
            weight.*exp(1i*fitError))/max(sum(weight),eps));
        edgeWeight(ip) = 1/max(fitRmse^2,deg2rad(1)^2);
    end

    % Independent pair fits may differ by integer multiples of 2*pi.  Align
    % them to the spanning star rooted at sensor 1 before the graph solve.
    edgePhase = alignPairBranches(edgePhase,pairIndex,nSensor);

    % Cap reliability contrast so one exceptionally clean edge cannot erase
    % the remaining sensor evidence.  For two sensors this reduces exactly
    % to the familiar half-phase, zero-mean differential gauge.
    edgeWeight = edgeWeight/max(median(edgeWeight),eps);
    edgeWeight = min(max(edgeWeight,0.1),10);
    graph.edgeWeight(:,ih) = edgeWeight;
    activeEdge = graph.edgeFitCoherence(:,ih)>=minEdgeFitCoherence;
    graph.activeEdgeMask(:,ih) = activeEdge;
    [relativePhase,activeSensor] = projectReliableComponents( ...
        edgePhase,pairIndex,incidence,edgeWeight,activeEdge,nSensor);
    graph.activeSensorMask(:,ih) = activeSensor;

    projectedEdge = relativePhase*incidence.';
    if any(activeEdge)
        projectionError = angle(exp(1i*( ...
            edgePhase(:,activeEdge)-projectedEdge(:,activeEdge))));
        activeWeight = edgeWeight(activeEdge);
        graph.edgeProjectionRmseDeg(ih) = rad2deg(sqrt( ...
            sum(projectionError.^2*activeWeight)/ ...
            max(nSample*sum(activeWeight),eps)));
    else
        graph.edgeProjectionRmseDeg(ih) = nan;
    end
    graph.cycleClosureRmseDeg(ih) = cycleClosureRmse( ...
        edgePhase,pairIndex,nSensor,activeEdge);
    for is = 1:nSensor
        pathPhase = relativePhase(:,is);
        pathPhasor(:,is,ih) = exp(1i*pathPhase);
        corrected(:,is,ih) = component(:,is,ih).*exp(-1i*pathPhase);
    end
end
graph.meanCycleClosureRmseDeg = mean( ...
    graph.cycleClosureRmseDeg,'omitnan');
graph.minEdgeFitCoherence = minEdgeFitCoherence;
graph.definition = [ ...
    'All pairwise fitted path phases are projected onto zero-mean sensor ', ...
    'phase potentials. cycleClosureRmseDeg measures the inconsistency of ', ...
    'the independently fitted edges before graph projection.'];
end

function [relativePhase,activeSensor] = projectReliableComponents( ...
    edgePhase,pairIndex,incidence,edgeWeight,activeEdge,nSensor)
relativePhase = zeros(size(edgePhase,1),nSensor);
activeSensor = false(nSensor,1);
adjacency = false(nSensor);
for ip = find(activeEdge(:)).'
    adjacency(pairIndex(ip,1),pairIndex(ip,2)) = true;
    adjacency(pairIndex(ip,2),pairIndex(ip,1)) = true;
end
visited = false(nSensor,1);
for startSensor = 1:nSensor
    if visited(startSensor)
        continue;
    end
    queue = startSensor;
    componentSensor = [];
    visited(startSensor) = true;
    while ~isempty(queue)
        sensor = queue(1);
        queue(1) = [];
        componentSensor(end+1) = sensor; %#ok<AGROW>
        neighbor = find(adjacency(sensor,:) & ~visited.');
        visited(neighbor) = true;
        queue = [queue,neighbor]; %#ok<AGROW>
    end
    if numel(componentSensor)<2
        continue;
    end
    componentEdge = find(activeEdge & ...
        ismember(pairIndex(:,1),componentSensor) & ...
        ismember(pairIndex(:,2),componentSensor));
    localIncidence = incidence(componentEdge,componentSensor);
    localWeight = edgeWeight(componentEdge);
    laplacian = localIncidence'*diag(localWeight)*localIncidence;
    graphMap = (laplacian+ones(numel(componentSensor))/ ...
        numel(componentSensor))\(localIncidence'*diag(localWeight));
    relativePhase(:,componentSensor) = ...
        edgePhase(:,componentEdge)*graphMap.';
    activeSensor(componentSensor) = true;
end
end

function [pairIndex,incidence] = completePairGraph(nSensor)
nPair = nSensor*(nSensor-1)/2;
pairIndex = zeros(nPair,2);
incidence = zeros(nPair,nSensor);
ip = 0;
for sensorI = 1:nSensor-1
    for sensorJ = sensorI+1:nSensor
        ip = ip+1;
        pairIndex(ip,:) = [sensorI,sensorJ];
        incidence(ip,sensorI) = -1;
        incidence(ip,sensorJ) = 1;
    end
end
end

function edgePhase = alignPairBranches(edgePhase,pairIndex,nSensor)
referencePhase = zeros(size(edgePhase,1),nSensor);
for sensorJ = 2:nSensor
    edge = find(pairIndex(:,1)==1 & pairIndex(:,2)==sensorJ,1);
    referencePhase(:,sensorJ) = edgePhase(:,edge);
end
for ip = 1:size(pairIndex,1)
    sensorI = pairIndex(ip,1);
    sensorJ = pairIndex(ip,2);
    referenceDifference = referencePhase(:,sensorJ)- ...
        referencePhase(:,sensorI);
    branchShift = 2*pi*round(median( ...
        (referenceDifference-edgePhase(:,ip))/(2*pi)));
    edgePhase(:,ip) = edgePhase(:,ip)+branchShift;
end
end

function value = cycleClosureRmse(edgePhase,pairIndex,nSensor,activeEdge)
closure = [];
for sensorI = 1:nSensor-2
    for sensorJ = sensorI+1:nSensor-1
        for sensorK = sensorJ+1:nSensor
            edgeIJ = find(pairIndex(:,1)==sensorI & ...
                pairIndex(:,2)==sensorJ,1);
            edgeJK = find(pairIndex(:,1)==sensorJ & ...
                pairIndex(:,2)==sensorK,1);
            edgeIK = find(pairIndex(:,1)==sensorI & ...
                pairIndex(:,2)==sensorK,1);
            if ~all(activeEdge([edgeIJ,edgeJK,edgeIK]))
                continue;
            end
            triangle = angle(exp(1i*(edgePhase(:,edgeIJ)+ ...
                edgePhase(:,edgeJK)-edgePhase(:,edgeIK))));
            closure = [closure;triangle]; %#ok<AGROW>
        end
    end
end
if isempty(closure)
    value = nan;
else
    value = rad2deg(sqrt(mean(closure.^2)));
end
end

function weight = reliabilityWeights(component)
amplitude = abs(component);
[~,nSensor,nHarmonic] = size(component);
weight = zeros(size(amplitude));
for is = 1:nSensor
    for ih = 1:nHarmonic
        scale = median(amplitude(:,is,ih));
        normalized = amplitude(:,is,ih)/max(scale,eps);
        weight(:,is,ih) = normalized.^2./(1+normalized.^2);
    end
end
end

function common = incrementalCommonPhase(component,harmonicOrder, ...
    reference,weight,fs)
[nSample,nSensor,nHarmonic] = size(component);
referenceIncrement = diff(reference);
numerator = zeros(nSample-1,1);
denominator = zeros(nSample-1,1);
for is = 1:nSensor
    for ih = 1:nHarmonic
        observedIncrement = angle(component(2:end,is,ih).* ...
            conj(component(1:end-1,is,ih)).* ...
            exp(-1i*harmonicOrder(ih)*referenceIncrement));
        observedIncrement = observedIncrement/harmonicOrder(ih);
        pairWeight = min(weight(2:end,is,ih),weight(1:end-1,is,ih));
        numerator = numerator+pairWeight.*observedIncrement;
        denominator = denominator+pairWeight;
    end
end
correctionIncrement = numerator./max(denominator,eps);
correctionIncrement = movingAverage(correctionIncrement, ...
    max(3,round(0.0015*fs)));
common = reference+[0;cumsum(correctionIncrement)];
end

function phase = regularizeCommonPhase(raw,nominal,carrierPhase,fs,options)
window = max(3,round(options.phaseSmoothSec*fs));
delta = movingAverage(raw-nominal,window);
if options.removeCommonSpatialOrder && options.pathOrder>0
    positiveOrder = options.spatialMultiple*(1:options.pathOrder);
    design = ones(numel(delta),1+2*numel(positiveOrder));
    for ik = 1:numel(positiveOrder)
        design(:,2*ik) = cos(positiveOrder(ik)*carrierPhase);
        design(:,2*ik+1) = sin(positiveOrder(ik)*carrierPhase);
    end
    coefficient = design\delta;
    % Retain the intercept but remove carrier-synchronous common phase.
    delta = delta-design(:,2:end)*coefficient(2:end);
end
phase = nominal+delta;
end

function coefficient = weightedComplexFit(basis,target,weight,ridge,penalty)
rootWeight = sqrt(max(weight,0));
weightedBasis = basis.*rootWeight;
weightedTarget = target.*rootWeight;
normal = weightedBasis'*weightedBasis;
coefficient = (normal+ridge*trace(normal)/size(normal,1)* ...
    diag(penalty))\(weightedBasis'*weightedTarget);
end

function y = movingAverage(x,window)
window = max(1,2*floor(window/2)+1);
kernel = ones(window,1)/window;
pad = floor(window/2);
xp = [repmat(x(1,:),pad,1);x;repmat(x(end,:),pad,1)];
y = conv2(xp,kernel,'valid');
end
