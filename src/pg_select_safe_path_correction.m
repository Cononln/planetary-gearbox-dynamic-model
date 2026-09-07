function result = pg_select_safe_path_correction(component,carrierPhase, ...
    spatialOrder,pathPenalty,options)
%PG_SELECT_SAFE_PATH_CORRECTION Cross-validated path-model safeguard.
%
% Candidate moving-path corrections are trained on alternating carrier
% revolutions and evaluated on the revolutions that were not used for the
% fit.  The candidates are no correction, each single sensor-pair edge,
% and the complete sensor graph.  A path model is accepted only when its
% held-out phase coherence improves by a prescribed margin without a bad
% fold.  The winning topology is then refitted using the full record.
%
% The returned correction is unit modulus and therefore cannot change the
% analytic-component amplitude.  If no candidate passes the safeguard,
% the input is returned unchanged.

if nargin < 5
    options = struct;
end
options = fillDefaults(options);
carrierPhase = carrierPhase(:);
spatialOrder = spatialOrder(:).';
pathPenalty = pathPenalty(:).';
[nSample,nSensor,nHarmonic] = size(component);
if numel(carrierPhase)~=nSample
    error('pg_select_safe_path_correction:PhaseSize', ...
        'carrierPhase must contain one value per component sample.');
end
if numel(spatialOrder)~=numel(pathPenalty)
    error('pg_select_safe_path_correction:PenaltySize', ...
        'spatialOrder and pathPenalty must have the same length.');
end
if nSensor<2
    error('pg_select_safe_path_correction:SensorCount', ...
        'At least two synchronized sensors are required.');
end

pairIndex = nchoosek(1:nSensor,2);
nPair = size(pairIndex,1);
[candidateName,candidatePairMask,candidateSensorMask] = ...
    makeCandidates(pairIndex,nSensor);
nCandidate = numel(candidateName);

stride = max(1,round(options.cvStride));
cvIndex = (1:stride:nSample).';
cvComponent = component(cvIndex,:,:);
cvCarrierPhase = carrierPhase(cvIndex);
cvBasis = exp(1i*cvCarrierPhase*spatialOrder);
carrierCycle = floor((cvCarrierPhase-cvCarrierPhase(1))/(2*pi));
foldId = mod(carrierCycle,2)+1;
if numel(unique(foldId))<2
    error('pg_select_safe_path_correction:CarrierCycles', ...
        'The record must span at least two carrier revolutions.');
end

rawCoherence = nan(nCandidate,2);
correctedCoherence = nan(nCandidate,2);
gain = zeros(nCandidate,2);
trainEdgeCoherence = nan(nCandidate,2);
for ic = 1:nCandidate
    activePair = candidatePairMask(ic,:).';
    if ~any(activePair)
        for iFold = 1:2
            test = foldId==iFold;
            rawCoherence(ic,iFold) = graphCoherence( ...
                cvComponent(test,:,:),pairIndex, ...
                true(nPair,nHarmonic),[]);
            correctedCoherence(ic,iFold) = rawCoherence(ic,iFold);
        end
        continue;
    end
    activePairHarmonic = repmat(activePair,1,nHarmonic);
    for iFold = 1:2
        test = foldId==iFold;
        train = ~test;
        [predictedEdge,edgeCoherence] = fitEdges( ...
            cvComponent,cvBasis,train,pairIndex,activePairHarmonic, ...
            options.ridge,pathPenalty);
        projectedEdge = projectEdges(predictedEdge,pairIndex, ...
            activePairHarmonic,edgeCoherence,nSensor);
        rawCoherence(ic,iFold) = graphCoherence( ...
            cvComponent(test,:,:),pairIndex, ...
            activePairHarmonic,[]);
        correctedCoherence(ic,iFold) = graphCoherence( ...
            cvComponent(test,:,:),pairIndex, ...
            activePairHarmonic,projectedEdge(test,:,:));
        gain(ic,iFold) = correctedCoherence(ic,iFold)- ...
            rawCoherence(ic,iFold);
        trainEdgeCoherence(ic,iFold) = mean( ...
            edgeCoherence(activePairHarmonic),'omitnan');
    end
end

meanGain = mean(gain,2);
worstFoldGain = min(gain,[],2);
meanCorrectedCoherence = mean(correctedCoherence,2);
gainFoldSpread = abs(gain(:,1)-gain(:,2));
correctedFoldSpread = abs(correctedCoherence(:,1)- ...
    correctedCoherence(:,2));
% Positive gain is an admission rule, not the ranking objective.  Ranking
% by gain alone favors a very poor edge that merely starts near zero.  Once
% admitted, choose the most coherent and fold-stable held-out topology.
selectionScore = meanCorrectedCoherence- ...
    options.stabilityPenalty*correctedFoldSpread;
admissible = meanGain>=options.minCvGain & ...
    worstFoldGain>=-options.maxFoldLoss;
admissible(1) = true;
selectionScore(~admissible) = -inf;
[~,selectedCandidate] = max(selectionScore);

selectedPairMask = candidatePairMask(selectedCandidate,:).';
selectedSensorMask = candidateSensorMask(selectedCandidate,:).';
if any(selectedPairMask)
    basis = exp(1i*carrierPhase*spatialOrder);
    activePairHarmonic = repmat(selectedPairMask,1,nHarmonic);
    [predictedEdge,edgeCoherence] = fitEdges( ...
        component,basis,true(nSample,1),pairIndex, ...
        activePairHarmonic,options.ridge,pathPenalty);
    projectedEdge = projectEdges(predictedEdge,pairIndex, ...
        activePairHarmonic,edgeCoherence,nSensor);
    [pathPhasor,correctedComponent] = applyProjectedEdges( ...
        component,projectedEdge,pairIndex,activePairHarmonic, ...
        edgeCoherence,nSensor,selectedSensorMask);
    fullRawCoherence = graphCoherence(component,pairIndex, ...
        activePairHarmonic,[]);
    fullCorrectedCoherence = graphCoherence(component,pairIndex, ...
        activePairHarmonic,projectedEdge);
    fullGain = fullCorrectedCoherence-fullRawCoherence;
else
    edgeCoherence = zeros(nPair,nHarmonic);
    pathPhasor = complex(ones(size(component)));
    correctedComponent = component;
    fullRawCoherence = graphCoherence(component,pairIndex, ...
        true(nPair,nHarmonic),[]);
    fullCorrectedCoherence = fullRawCoherence;
    fullGain = 0;
end

% This final guard protects against numerical or full-record refit failure.
if fullGain < -options.fullFitTolerance
    selectedCandidate = 1;
    selectedPairMask = false(nPair,1);
    selectedSensorMask = true(nSensor,1);
    edgeCoherence = zeros(nPair,nHarmonic);
    pathPhasor = complex(ones(size(component)));
    correctedComponent = component;
    fullCorrectedCoherence = graphCoherence(component,pairIndex, ...
        true(nPair,nHarmonic),[]);
    fullRawCoherence = fullCorrectedCoherence;
    fullGain = 0;
end

% Optional cross-fitted output for an independent downstream validation.
% Every sample is corrected by a model trained on the opposite carrier-
% revolution parity.  The default is off because full-rate measured data
% can be large; validation runners enable it after band-limited resampling.
selectedSensorIndex = find(selectedSensorMask).';
crossFittedCorrectedSelected = [];
crossFitFoldId = [];
if options.returnCrossFittedComponent
    fullCarrierCycle = floor((carrierPhase-carrierPhase(1))/(2*pi));
    crossFitFoldId = mod(fullCarrierCycle,2)+1;
    if any(selectedPairMask)
        basis = exp(1i*carrierPhase*spatialOrder);
        activePairHarmonic = repmat(selectedPairMask,1,nHarmonic);
        crossFittedCorrectedSelected = complex(zeros( ...
            nSample,numel(selectedSensorIndex),nHarmonic));
        for iFold = 1:2
            test = crossFitFoldId==iFold;
            train = ~test;
            [foldPredictedEdge,foldEdgeCoherence] = fitEdges( ...
                component,basis,train,pairIndex,activePairHarmonic, ...
                options.ridge,pathPenalty);
            foldProjectedEdge = projectEdges(foldPredictedEdge, ...
                pairIndex,activePairHarmonic,foldEdgeCoherence,nSensor);
            [~,foldCorrected] = applyProjectedEdges(component, ...
                foldProjectedEdge,pairIndex,activePairHarmonic, ...
                foldEdgeCoherence,nSensor,selectedSensorMask);
            crossFittedCorrectedSelected(test,:,:) = ...
                foldCorrected(test,selectedSensorIndex,:);
        end
    else
        crossFittedCorrectedSelected = ...
            component(:,selectedSensorIndex,:);
    end
end

sensorSubset = strings(nCandidate,1);
edgeSubset = strings(nCandidate,1);
for ic = 1:nCandidate
    sensorSubset(ic) = strjoin(string(find( ...
        candidateSensorMask(ic,:))),'+');
    active = find(candidatePairMask(ic,:));
    if isempty(active)
        edgeSubset(ic) = "none";
    else
        label = strings(numel(active),1);
        for k = 1:numel(active)
            label(k) = sprintf('%d-%d',pairIndex(active(k),1), ...
                pairIndex(active(k),2));
        end
        edgeSubset(ic) = strjoin(label,'+');
    end
end
selected = false(nCandidate,1);
selected(selectedCandidate) = true;
cvTable = table(candidateName,sensorSubset,edgeSubset, ...
    rawCoherence(:,1),correctedCoherence(:,1),gain(:,1), ...
    rawCoherence(:,2),correctedCoherence(:,2),gain(:,2), ...
    mean(trainEdgeCoherence,2,'omitnan'),meanGain,worstFoldGain, ...
    meanCorrectedCoherence,gainFoldSpread,correctedFoldSpread, ...
    selectionScore,admissible,selected, ...
    'VariableNames',{'candidate','sensor_subset','edge_subset', ...
    'fold1_raw_coherence','fold1_corrected_coherence','fold1_gain', ...
    'fold2_raw_coherence','fold2_corrected_coherence','fold2_gain', ...
    'mean_train_edge_coherence','mean_gain','worst_fold_gain', ...
    'mean_corrected_coherence','gain_fold_spread', ...
    'corrected_fold_spread','selection_score','admissible','selected'});

result.selectedCandidate = candidateName(selectedCandidate);
result.selectedCandidateIndex = selectedCandidate;
result.selectedSensorMask = selectedSensorMask;
result.selectedSensorIndex = selectedSensorIndex;
result.selectedPairMask = selectedPairMask;
result.pathPhasor = pathPhasor;
result.correctedComponent = correctedComponent;
result.correctedSelectedComponent = ...
    correctedComponent(:,selectedSensorIndex,:);
result.selectedPathPhasor = pathPhasor(:,selectedSensorIndex,:);
result.crossFittedCorrectedSelectedComponent = ...
    crossFittedCorrectedSelected;
result.crossFitFoldId = crossFitFoldId;
result.cvTable = cvTable;
result.fullRawCoherence = fullRawCoherence;
result.fullCorrectedCoherence = fullCorrectedCoherence;
result.fullCoherenceGain = fullGain;
result.amplitudePreservationError = max(abs(abs(correctedComponent(:))- ...
    abs(component(:))))/max(abs(component(:)));
result.options = options;
result.pathGraph = makeSelectedGraph(pairIndex,selectedPairMask, ...
    selectedSensorIndex,edgeCoherence,nHarmonic);
result.decisionRule = [ ...
    'Alternating carrier-revolution two-fold validation. Accept a path ', ...
    'topology only when mean held-out coherence gain exceeds minCvGain ', ...
    'and neither fold loses more than maxFoldLoss. Among admissible ', ...
    'topologies choose the highest fold-stable corrected coherence; ', ...
    'otherwise return the unmodified components.'];
end

function options = fillDefaults(options)
defaults.cvStride = 20;
defaults.ridge = 2e-3;
defaults.minCvGain = 0.015;
defaults.maxFoldLoss = 0.010;
defaults.stabilityPenalty = 0.25;
defaults.fullFitTolerance = 1e-6;
defaults.returnCrossFittedComponent = false;
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(options,names{i}) || isempty(options.(names{i}))
        options.(names{i}) = defaults.(names{i});
    end
end
end

function [name,pairMask,sensorMask] = makeCandidates(pairIndex,nSensor)
nPair = size(pairIndex,1);
nCandidate = nPair+2;
name = strings(nCandidate,1);
pairMask = false(nCandidate,nPair);
sensorMask = false(nCandidate,nSensor);
name(1) = "no_path";
sensorMask(1,:) = true;
for ip = 1:nPair
    name(ip+1) = sprintf('pair_%d_%d',pairIndex(ip,1), ...
        pairIndex(ip,2));
    pairMask(ip+1,ip) = true;
    sensorMask(ip+1,pairIndex(ip,:)) = true;
end
name(end) = "complete_graph";
pairMask(end,:) = true;
sensorMask(end,:) = true;
end

function [predictedEdge,fitCoherence] = fitEdges( ...
    component,basis,train,pairIndex,activePairHarmonic,ridge,pathPenalty)
[nSample,~,nHarmonic] = size(component);
nPair = size(pairIndex,1);
predictedEdge = complex(ones(nSample,nPair,nHarmonic));
fitCoherence = zeros(nPair,nHarmonic);
for ih = 1:nHarmonic
    for ip = 1:nPair
        if ~activePairHarmonic(ip,ih)
            continue;
        end
        sensorI = pairIndex(ip,1);
        sensorJ = pairIndex(ip,2);
        cross = component(:,sensorJ,ih).*conj(component(:,sensorI,ih));
        unitCross = cross./max(abs(cross),eps);
        amplitude = sqrt(abs(component(:,sensorI,ih)).* ...
            abs(component(:,sensorJ,ih)));
        normalized = amplitude/max(median(amplitude(train)),eps);
        weight = normalized.^2./(1+normalized.^2);
        coefficient = complexFit(basis(train,:),unitCross(train), ...
            weight(train),ridge,pathPenalty);
        prediction = basis*coefficient;
        prediction = prediction./max(abs(prediction),eps);
        predictedEdge(:,ip,ih) = prediction;
        error = angle(unitCross(train).*conj(prediction(train)));
        fitCoherence(ip,ih) = abs(sum(weight(train).*exp(1i*error))/ ...
            max(sum(weight(train)),eps));
    end
end
end

function projectedEdge = projectEdges(predictedEdge,pairIndex, ...
    activePairHarmonic,edgeCoherence,nSensor)
[nSample,nPair,nHarmonic] = size(predictedEdge);
projectedEdge = complex(ones(nSample,nPair,nHarmonic));
incidence = zeros(nPair,nSensor);
for ip = 1:nPair
    incidence(ip,pairIndex(ip,1)) = -1;
    incidence(ip,pairIndex(ip,2)) = 1;
end
for ih = 1:nHarmonic
    active = activePairHarmonic(:,ih);
    if ~any(active)
        continue;
    end
    edgePhase = zeros(nSample,nPair);
    edgePhase(:,active) = unwrap(angle(predictedEdge(:,active,ih)));
    edgePhase = alignActiveBranches(edgePhase,pairIndex,active,nSensor);
    activeSensor = unique(pairIndex(active,:));
    localIncidence = incidence(active,activeSensor);
    weight = max(edgeCoherence(active,ih).^2,0.05);
    laplacian = localIncidence'*diag(weight)*localIncidence;
    graphMap = (laplacian+ones(numel(activeSensor))/ ...
        numel(activeSensor))\(localIncidence'*diag(weight));
    sensorPhase = edgePhase(:,active)*graphMap.';
    projected = sensorPhase*localIncidence.';
    projectedEdge(:,active,ih) = exp(1i*projected);
end
end

function edgePhase = alignActiveBranches(edgePhase,pairIndex,active,nSensor)
activeIndex = find(active);
if numel(activeIndex)<2
    return;
end
reference = nan(size(edgePhase,1),nSensor);
reference(:,1) = 0;
updated = true;
while updated
    updated = false;
    for ip = activeIndex(:).'
        sensorI = pairIndex(ip,1);
        sensorJ = pairIndex(ip,2);
        if all(~isnan(reference(:,sensorI))) && ...
                any(isnan(reference(:,sensorJ)))
            reference(:,sensorJ) = reference(:,sensorI)+edgePhase(:,ip);
            updated = true;
        elseif all(~isnan(reference(:,sensorJ))) && ...
                any(isnan(reference(:,sensorI)))
            reference(:,sensorI) = reference(:,sensorJ)-edgePhase(:,ip);
            updated = true;
        end
    end
end
for ip = activeIndex(:).'
    sensorI = pairIndex(ip,1);
    sensorJ = pairIndex(ip,2);
    if all(~isnan(reference(:,sensorI))) && ...
            all(~isnan(reference(:,sensorJ)))
        target = reference(:,sensorJ)-reference(:,sensorI);
        edgePhase(:,ip) = edgePhase(:,ip)+2*pi*round( ...
            median((target-edgePhase(:,ip))/(2*pi)));
    end
end
end

function [pathPhasor,corrected] = applyProjectedEdges(component, ...
    projectedEdge,pairIndex,activePairHarmonic,edgeCoherence,nSensor, ...
    selectedSensorMask)
[nSample,~,nHarmonic] = size(component);
nPair = size(pairIndex,1);
pathPhasor = complex(ones(size(component)));
corrected = component;
incidence = zeros(nPair,nSensor);
for ip = 1:nPair
    incidence(ip,pairIndex(ip,1)) = -1;
    incidence(ip,pairIndex(ip,2)) = 1;
end
for ih = 1:nHarmonic
    active = activePairHarmonic(:,ih);
    if ~any(active)
        continue;
    end
    activeSensor = find(selectedSensorMask).';
    localIncidence = incidence(active,activeSensor);
    edgePhase = unwrap(angle(projectedEdge(:,active,ih)));
    weight = max(edgeCoherence(active,ih).^2,0.05);
    laplacian = localIncidence'*diag(weight)*localIncidence;
    graphMap = (laplacian+ones(numel(activeSensor))/ ...
        numel(activeSensor))\(localIncidence'*diag(weight));
    sensorPhase = edgePhase*graphMap.';
    for k = 1:numel(activeSensor)
        sensor = activeSensor(k);
        pathPhasor(:,sensor,ih) = exp(1i*sensorPhase(:,k));
        corrected(:,sensor,ih) = component(:,sensor,ih).* ...
            conj(pathPhasor(:,sensor,ih));
    end
end
end

function value = graphCoherence(component,pairIndex,activePairHarmonic, ...
    predictedEdge)
[~,~,nHarmonic] = size(component);
sumValue = 0;
sumWeight = 0;
for ih = 1:nHarmonic
    for ip = 1:size(pairIndex,1)
        if ~activePairHarmonic(ip,ih)
            continue;
        end
        sensorI = pairIndex(ip,1);
        sensorJ = pairIndex(ip,2);
        cross = component(:,sensorJ,ih).*conj(component(:,sensorI,ih));
        unitCross = cross./max(abs(cross),eps);
        if ~isempty(predictedEdge)
            unitCross = unitCross.*conj(predictedEdge(:,ip,ih));
        end
        amplitude = sqrt(abs(component(:,sensorI,ih)).* ...
            abs(component(:,sensorJ,ih)));
        normalized = amplitude/max(median(amplitude),eps);
        weight = normalized.^2./(1+normalized.^2);
        coherence = abs(sum(weight.*unitCross)/max(sum(weight),eps));
        pairWeight = sqrt(sum(weight));
        sumValue = sumValue+pairWeight*coherence;
        sumWeight = sumWeight+pairWeight;
    end
end
value = sumValue/max(sumWeight,eps);
end

function coefficient = complexFit(basis,target,weight,ridge,penalty)
rootWeight = sqrt(max(weight,0));
weightedBasis = basis.*rootWeight;
normal = weightedBasis'*weightedBasis;
scale = max(real(trace(normal))/size(normal,1),eps);
coefficient = (normal+ridge*scale*diag(penalty))\ ...
    (weightedBasis'*(target.*rootWeight));
end

function graph = makeSelectedGraph(globalPairIndex,selectedPairMask, ...
    selectedSensorIndex,edgeCoherence,nHarmonic)
nSelected = numel(selectedSensorIndex);
if ~any(selectedPairMask)
    graph = [];
    return;
end
graph.pairIndex = nchoosek(1:nSelected,2);
graph.activeEdgeMask = true(size(graph.pairIndex,1),nHarmonic);
graph.activeSensorMask = true(nSelected,nHarmonic);
graph.edgeFitCoherence = zeros(size(graph.pairIndex,1),nHarmonic);
for ip = 1:size(graph.pairIndex,1)
    globalPair = sort(selectedSensorIndex(graph.pairIndex(ip,:)));
    index = find(globalPairIndex(:,1)==globalPair(1) & ...
        globalPairIndex(:,2)==globalPair(2),1);
    graph.edgeFitCoherence(ip,:) = edgeCoherence(index,:);
end
graph.cycleClosureRmseDeg = nan(nHarmonic,1);
graph.meanCycleClosureRmseDeg = nan;
graph.definition = ['Topology selected by alternating-carrier-cycle ', ...
    'held-out coherence; local pair indices refer to selectedSensorIndex.'];
end
