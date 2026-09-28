function [H,energy,detail] = pg_fe_ring_modal_path_gain(theta,frequencyHz,model,options)
%PG_FE_RING_MODAL_PATH_GAIN Modal-shape-only contact-to-sensor transfer.
%
%   H(theta,sensor,frequency) is the complex radial acceleration at an
%   exported sensor per unit normal force applied at the moving ring
%   contact angle theta.  The transfer is assembled directly from the
%   exported ANSYS input and sensor mode-shape coefficients:
%
%       H_s(theta,w) = sum_m psi_in(theta,m) *
%           [-w^2/(wn_m^2-w^2+i*2*zeta_m*wn_m*w)] * psi_s(s,m).
%
%   No angular-distance gain, fitted sensor gain, or empirical attenuation
%   is used.  The optional ``energy`` output is |H|^2.  It is useful for
%   diagnosing modal participation/constraint attenuation, not for
%   replacing the time-domain path solver.
%
%   theta may be a scalar or vector (rad), frequencyHz may be a scalar or
%   vector (Hz).  H has size [numel(theta), nSensor, numel(frequencyHz)].

arguments
    theta (:,1) double
    frequencyHz (:,1) double
    model struct
    options.DampingRatio (:,1) double {mustBePositive} = model.zeta(:)
end
theta = mod(theta(:),2*pi);
frequencyHz = frequencyHz(:);
nTheta = numel(theta);
nSensor = size(model.sensorShapeRadial,1);
nMode = numel(model.frequencyHz);
if size(model.inputShapeNormal,2) ~= nMode || ...
        size(model.sensorShapeRadial,2) ~= nMode
    error('pg_fe_ring_modal_path_gain:ModeSize', ...
        'Input and sensor shape mode counts are inconsistent.');
end
zeta = options.DampingRatio(:).';
if numel(zeta) ~= nMode
    error('pg_fe_ring_modal_path_gain:DampingSize', ...
        'DampingRatio must have one value per retained mode.');
end

% Periodic interpolation of the exported moving-contact participation.
a = mod(model.contactAngleRad(:),2*pi);
[a,order] = sort(a);
psi = model.inputShapeNormal(order,:);
psiTheta = zeros(nTheta,nMode);
for im = 1:nMode
    psiTheta(:,im) = interp1([a(end)-2*pi;a;a(1)+2*pi], ...
        [psi(end,im);psi(:,im);psi(1,im)],theta,'linear');
end

wn = 2*pi*model.frequencyHz(:).';
H = complex(zeros(nTheta,nSensor,numel(frequencyHz)));
for jf = 1:numel(frequencyHz)
    w = 2*pi*frequencyHz(jf);
    modalH = -(w.^2)./(wn.^2-w.^2 + 1i*2*zeta.*wn*w);
    % A contact-by-sensor matrix for this frequency.
    H(:,:,jf) = (psiTheta.*modalH) * model.sensorShapeRadial.';
end
energy = abs(H).^2;
detail = struct();
detail.thetaRad = theta;
detail.frequencyHz = frequencyHz;
detail.sensorAnglesLocalDeg = model.sensorAnglesLocalDeg(:).';
detail.modeFrequencyHz = model.frequencyHz(:).';
detail.dampingRatio = zeta;
detail.inputParticipation = psiTheta;
detail.sensorShapeRadial = model.sensorShapeRadial;
detail.description = ['Modal-shape-only transfer; no empirical distance ' ...
    'weight.  Energy is |H|^2 per unit squared normal force.'];
end
