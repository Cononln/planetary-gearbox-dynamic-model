function model = pg_ring_modal_model(p)
%PG_RING_MODAL_MODEL Analytical ring model calibrated by supplied frequencies.
% Generalized coordinates eta have units of radial displacement (m).

if nargin < 1
    p = pg_parameters;
end
thisFile = mfilename('fullpath');
projectRoot = fileparts(fileparts(thisFile));
modeTable = readtable(fullfile(projectRoot,'data','ring_modes_30.csv'));
rows = ismember(modeTable.mode,p.ringModal.sourceModes);
frequencyHz = modeTable.frequency_hz(rows);
sourceModes = modeTable.mode(rows);
[sourceModes,order] = sort(sourceModes);
frequencyHz = frequencyHz(order);

numberOfModes = numel(p.ringModal.sourceModes);
if numel(frequencyHz) ~= numberOfModes
    error('pg_ring_modal_model:ModeCount', ...
        'Expected %d supplied ring modes.',numberOfModes);
end
if numel(p.ringModal.waveNumbers) ~= numberOfModes || ...
        numel(p.ringModal.family) ~= numberOfModes
    error('pg_ring_modal_model:ModeDefinition', ...
        'sourceModes, waveNumbers, and family must have equal lengths.');
end

model.source = 'frequency-calibrated analytical ring basis';
model.sourceModes = sourceModes(:).';
model.frequencyHz = frequencyHz(:).';
model.omega = 2*pi*model.frequencyHz;
model.waveNumber = p.ringModal.waveNumbers;
model.family = p.ringModal.family;
model.zeta = p.ringModal.zeta*ones(1,numberOfModes);

% With the displacement-normalized analytical shapes used here,
% circumference-mean psi^2 is 1/2 for n>0 sine/cosine modes and 1 for the
% n=0 cosine (uniform radial breathing) mode.
shapeMeanSquare = 0.5*ones(1,numberOfModes);
isAxisymmetric = model.waveNumber == 0;
if any(isAxisymmetric & ~strcmp(model.family,'cos'))
    error('pg_ring_modal_model:AxisymmetricFamily', ...
        'The n=0 breathing mode must use the cosine family.');
end
shapeMeanSquare(isAxisymmetric) = 1;
model.modalMass = p.body.ring.mass*shapeMeanSquare;
if isfield(p.ringModal,'modalMassScale')
    if numel(p.ringModal.modalMassScale) ~= numberOfModes || ...
            any(p.ringModal.modalMassScale <= 0)
        error('pg_ring_modal_model:ModalMassScale', ...
            'modalMassScale must contain one positive value per mode.');
    end
    model.modalMassScale = p.ringModal.modalMassScale;
else
    model.modalMassScale = ones(1,numberOfModes);
end
model.modalMass = model.modalMass.*model.modalMassScale;
model.M = diag(model.modalMass);
model.C = diag(2*model.zeta.*model.omega.*model.modalMass);
model.K = diag((model.omega.^2).*model.modalMass);
model.radialMeshProjection = p.ringModal.radialMeshProjection;
end
