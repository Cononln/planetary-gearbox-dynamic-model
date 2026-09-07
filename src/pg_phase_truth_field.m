function truth = pg_phase_truth_field(frequencyHz,phi_c,p,model,sourceIndex)
%PG_PHASE_TRUTH_FIELD Compute dual-sensor phase and group-delay truth fields.
% Rows are frequency; columns are carrier angle. The source is one specified
% moving ring-planet mesh. Other planet paths are carrier-angle shifts.

if nargin < 5
    sourceIndex = 1;
end
frequencyHz = frequencyHz(:);
phi_c = phi_c(:).';
nF = numel(frequencyHz);
nPhi = numel(phi_c);
H0 = zeros(nF,nPhi);
H90 = zeros(nF,nPhi);
Cacc = pg_sensor_observation(p,model);
meanStiffness.sunPlanet = ...
    p.mesh.sunPlanet.kMean*ones(1,p.model.nPlanet);
meanStiffness.ringPlanet = ...
    p.mesh.ringPlanet.kMean*ones(1,p.model.nPlanet);

for iphi = 1:nPhi
    [M,C,K,meta] = pg_coupled_matrices(phi_c(iphi),p,model,meanStiffness);
    inputVector = meta.coupling(sourceIndex).g;
    for jf = 1:nF
        omega = 2*pi*frequencyHz(jf);
        D = K-omega^2*M+1i*omega*C;
        response = -omega^2*Cacc*(D\inputVector);
        H0(jf,iphi) = response(1);
        H90(jf,iphi) = response(2);
    end
end

phaseWrapped = angle(H90.*conj(H0));
phaseUnwrappedFrequency = unwrap(phaseWrapped,[],1);
omegaVector = 2*pi*frequencyHz;
groupDelay = zeros(size(phaseUnwrappedFrequency));
for iphi = 1:nPhi
    groupDelay(:,iphi) = -gradient(phaseUnwrappedFrequency(:,iphi),omegaVector);
end

amp0 = abs(H0);
amp90 = abs(H90);
balance = 2*amp0.*amp90./(amp0.^2+amp90.^2+eps);
relative0 = amp0./max(amp0,[],2);
relative90 = amp90./max(amp90,[],2);
reliability = balance.*sqrt(relative0.*relative90);
validMask = reliability>0.05 & relative0>1e-4 & relative90>1e-4;

bandsHz = [1750,1920;2920,3140;4890,5060];
bandDelay = nan(size(bandsHz,1),nPhi);
for ib = 1:size(bandsHz,1)
    inBand = frequencyHz>=bandsHz(ib,1) & frequencyHz<=bandsHz(ib,2);
    for iphi = 1:nPhi
        use = inBand & validMask(:,iphi) & isfinite(groupDelay(:,iphi));
        if nnz(use)>=3
            weight = reliability(use,iphi).*sqrt(amp0(use,iphi).*amp90(use,iphi));
            bandDelay(ib,iphi) = sum(weight.*groupDelay(use,iphi))/sum(weight);
        end
    end
end

truth.frequencyHz = frequencyHz;
truth.omega = omegaVector;
truth.phi_c = phi_c;
truth.sourceIndex = sourceIndex;
truth.H0 = H0;
truth.H90 = H90;
truth.phaseDifferenceWrapped = phaseWrapped;
truth.phaseDifferenceUnwrappedFrequency = phaseUnwrappedFrequency;
truth.groupDelaySec = groupDelay;
truth.reliability = reliability;
truth.validMask = validMask;
truth.bandsHz = bandsHz;
truth.bandDelaySec = bandDelay;
truth.definition = ['Delta phase = arg(H90)-arg(H0); ', ...
    'group delay = -d(Delta phase)/d(omega)'];
end
