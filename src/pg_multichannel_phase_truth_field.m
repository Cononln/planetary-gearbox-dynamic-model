function truth = pg_multichannel_phase_truth_field(frequencyHz,phi_c,p,model,sourceIndex,referenceSensor)
%PG_MULTICHANNEL_PHASE_TRUTH_FIELD Complex path truth for arbitrary sensors.
% H(f,phi,j) is the acceleration response at fixed sensor j from one
% moving ring-planet mesh source.  Differential phase and group delay are
% always reported relative to referenceSensor.

if nargin < 5 || isempty(sourceIndex)
    sourceIndex = 1;
end
sensor = pg_sensor_configuration(p);
if nargin < 6 || isempty(referenceSensor)
    referenceSensor = 1;
end
if ~isscalar(sourceIndex) || sourceIndex<1 || sourceIndex>p.model.nPlanet
    error('pg_multichannel_phase_truth_field:Source', ...
        'sourceIndex must identify one planet.');
end
if ~isscalar(referenceSensor) || referenceSensor<1 || ...
        referenceSensor>sensor.nSensor
    error('pg_multichannel_phase_truth_field:Reference', ...
        'referenceSensor must identify one sensor.');
end

frequencyHz = frequencyHz(:);
phi_c = phi_c(:).';
nFrequency = numel(frequencyHz);
nAngle = numel(phi_c);
nSensor = sensor.nSensor;
Cacc = pg_sensor_observation(p,model);
meanStiffness.sunPlanet = ...
    p.mesh.sunPlanet.kMean*ones(1,p.model.nPlanet);
meanStiffness.ringPlanet = ...
    p.mesh.ringPlanet.kMean*ones(1,p.model.nPlanet);

H = complex(zeros(nFrequency,nAngle,nSensor));
for ia = 1:nAngle
    [M,C,K,meta] = pg_coupled_matrices(phi_c(ia),p,model,meanStiffness);
    inputVector = meta.coupling(sourceIndex).g;
    for jf = 1:nFrequency
        omega = 2*pi*frequencyHz(jf);
        dynamicStiffness = K-omega^2*M+1i*omega*C;
        response = -omega^2*Cacc*(dynamicStiffness\inputVector);
        H(jf,ia,:) = reshape(response,[1,1,nSensor]);
    end
end

omega = 2*pi*frequencyHz;
phaseDifferenceWrapped = zeros(nFrequency,nAngle,nSensor);
phaseDifferenceUnwrapped = zeros(nFrequency,nAngle,nSensor);
groupDelaySec = zeros(nFrequency,nAngle,nSensor);
reliability = zeros(nFrequency,nAngle,nSensor);
validMask = false(nFrequency,nAngle,nSensor);
amplitude = abs(H);
amplitudeReference = amplitude(:,:,referenceSensor);
relativeReference = amplitudeReference./max(amplitudeReference,[],2);
for is = 1:nSensor
    phaseDifferenceWrapped(:,:,is) = angle(H(:,:,is).* ...
        conj(H(:,:,referenceSensor)));
    phaseDifferenceUnwrapped(:,:,is) = unwrap( ...
        phaseDifferenceWrapped(:,:,is),[],1);
    for ia = 1:nAngle
        groupDelaySec(:,ia,is) = -gradient( ...
            phaseDifferenceUnwrapped(:,ia,is),omega);
    end

    amplitudeTarget = amplitude(:,:,is);
    relativeTarget = amplitudeTarget./max(amplitudeTarget,[],2);
    balance = 2*amplitudeReference.*amplitudeTarget./ ...
        (amplitudeReference.^2+amplitudeTarget.^2+eps);
    reliability(:,:,is) = balance.*sqrt(relativeReference.*relativeTarget);
    validMask(:,:,is) = reliability(:,:,is)>0.05 & ...
        relativeReference>1e-4 & relativeTarget>1e-4;
end

truth.frequencyHz = frequencyHz;
truth.omega = omega;
truth.phi_c = phi_c;
truth.H = H;
truth.amplitude = amplitude;
truth.referenceSensor = referenceSensor;
truth.sourceIndex = sourceIndex;
truth.sensorAnglesDeg = sensor.locationAnglesDeg;
truth.sensorSensitiveAnglesDeg = sensor.sensitiveAnglesDeg;
truth.phaseDifferenceWrapped = phaseDifferenceWrapped;
truth.phaseDifferenceUnwrappedFrequency = phaseDifferenceUnwrapped;
truth.groupDelaySec = groupDelaySec;
truth.reliability = reliability;
truth.validMask = validMask;
truth.definition = [ ...
    'Differential phase and group delay are sensor j relative to the ' ...
    'reference sensor for a single moving ring-planet mesh source.'];
end
