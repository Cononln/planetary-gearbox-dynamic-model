function field = pg_path_field(frequencyHz,phi_c,p,model,sourceWeights)
%PG_PATH_FIELD Complex moving-mesh paths for all planets and both sensors.
% H(f,phi,sensor,planet) is the acceleration FRF from a unit equivalent
% ring-planet mesh force. HCombined is the coherent sum specified by the
% complex sourceWeights. Healthy in-phase meshes use [1 1 1].

if nargin < 5 || isempty(sourceWeights)
    sourceWeights = ones(1,p.model.nPlanet);
end
frequencyHz = frequencyHz(:);
phi_c = phi_c(:).';
sourceWeights = sourceWeights(:).';
if numel(sourceWeights) ~= p.model.nPlanet
    error('pg_path_field:SourceWeights', ...
        'sourceWeights must contain one complex weight per planet.');
end

nFrequency = numel(frequencyHz);
nAngle = numel(phi_c);
sensor = pg_sensor_configuration(p);
nSensor = sensor.nSensor;
nPlanet = p.model.nPlanet;
H = complex(zeros(nFrequency,nAngle,nSensor,nPlanet));
radialProjection = zeros(nAngle,nSensor,nPlanet);
tangentialProjection = zeros(nAngle,nSensor,nPlanet);

meanStiffness.sunPlanet = ...
    p.mesh.sunPlanet.kMean*ones(1,nPlanet);
meanStiffness.ringPlanet = ...
    p.mesh.ringPlanet.kMean*ones(1,nPlanet);
Cacc = pg_sensor_observation(p,model);

for ia = 1:nAngle
    fixedGeometry = pg_fixed_path_kinematics(phi_c(ia),p);
    radialProjection(ia,:,:) = reshape(fixedGeometry.radialProjection, ...
        [1,nSensor,nPlanet]);
    tangentialProjection(ia,:,:) = reshape( ...
        fixedGeometry.tangentialProjection,[1,nSensor,nPlanet]);
    [M,C,K,meta] = pg_coupled_matrices(phi_c(ia),p,model,meanStiffness);
    inputMatrix = zeros(size(M,1),nPlanet);
    for ip = 1:nPlanet
        inputMatrix(:,ip) = meta.coupling(ip).g;
    end
    for jf = 1:nFrequency
        omega = 2*pi*frequencyHz(jf);
        dynamicStiffness = K-omega^2*M+1i*omega*C;
        response = -omega^2*Cacc*(dynamicStiffness\inputMatrix);
        H(jf,ia,:,:) = reshape(response,[1,1,nSensor,nPlanet]);
    end
end

weights = reshape(sourceWeights,[1,1,1,nPlanet]);
HCombined = sum(H.*weights,4);
amplitude = abs(HCombined);
if nSensor >= 2
    % Backward-compatible differential phase: channel 2 relative to 1.
    phaseDifference = angle(HCombined(:,:,2).*conj(HCombined(:,:,1)));
    amplitudeReference = amplitude(:,:,1);
    amplitudeTarget = amplitude(:,:,2);
    relativeReference = amplitudeReference./max(amplitudeReference,[],2);
    relativeTarget = amplitudeTarget./max(amplitudeTarget,[],2);
    balance = 2*amplitudeReference.*amplitudeTarget./ ...
        (amplitudeReference.^2+amplitudeTarget.^2+eps);
    reliability = balance.*sqrt(relativeReference.*relativeTarget);
    validMask = reliability>0.05 & relativeReference>1e-4 & ...
        relativeTarget>1e-4;
else
    phaseDifference = zeros(nFrequency,nAngle);
    reliability = ones(nFrequency,nAngle);
    validMask = true(nFrequency,nAngle);
end

% Pairwise differential phases support three or more channels while the
% legacy phaseDifference field remains channel 2 minus channel 1.
phaseDifferenceToReference = zeros(nFrequency,nAngle,nSensor);
for is = 1:nSensor
    phaseDifferenceToReference(:,:,is) = angle( ...
        HCombined(:,:,is).*conj(HCombined(:,:,1)));
end

field.frequencyHz = frequencyHz;
field.phi_c = phi_c;
field.sensorAnglesDeg = sensor.locationAnglesDeg;
field.sensorSensitiveAnglesDeg = sensor.sensitiveAnglesDeg;
field.sourceWeights = sourceWeights;
field.H = H;
field.HCombined = HCombined;
field.fixedGeometry.radialProjection = radialProjection;
field.fixedGeometry.tangentialProjection = tangentialProjection;
field.amplitude = amplitude;
field.phaseDifference = phaseDifference;
field.phaseDifferenceToReference = phaseDifferenceToReference;
field.reliability = reliability;
field.validMask = validMask;
field.definition = [ ...
    'H(f,phi,sensor,planet): acceleration per unit equivalent mesh force; ' ...
    'HCombined is the coherent source-weighted planet sum.'];
end
