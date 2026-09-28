function model = pg_fe_ring_modal_model(assetPath,options)
%PG_FE_RING_MODAL_MODEL Load an ANSYS modal asset with a mass-normalization assumption.

arguments
    assetPath (1,1) string
    options.ModalDampingRatio (1,1) double {mustBePositive} = 0.02
end
if ~isfile(assetPath)
    error('pg_fe_ring_modal_model:MissingAsset','Missing asset: %s',assetPath);
end
a = load(assetPath);
required = {'frequency_Hz','contact_angle_rad','input_shape_normal', ...
    'sensor_angle_local_deg','sensor_angle_global_deg', ...
    'sensor_shape_radial','mass_normalized'};
for k = 1:numel(required)
    if ~isfield(a,required{k})
        error('pg_fe_ring_modal_model:MissingField', ...
            'Asset is missing %s.',required{k});
    end
end
if ~logical(a.mass_normalized(1))
    error('pg_fe_ring_modal_model:Normalization', ...
        'Only mass-normalized modes are supported.');
end
model.source = 'ANSYS 2023 R1 mass-normalized FE modes';
model.assetPath = char(assetPath);
model.frequencyHz = a.frequency_Hz(:).';
model.omega = 2*pi*model.frequencyHz;
model.zeta = options.ModalDampingRatio*ones(size(model.frequencyHz));
model.contactAngleRad = a.contact_angle_rad(:);
model.inputShapeNormal = a.input_shape_normal;
model.sensorAnglesLocalDeg = a.sensor_angle_local_deg(:).';
model.sensorAnglesGlobalDeg = a.sensor_angle_global_deg(:).';
model.sensorShapeRadial = a.sensor_shape_radial;
model.massNormalized = true;
model.massNormalizationSIVerified = isfield(a,'mass_normalization_si_verified') && ...
    logical(a.mass_normalization_si_verified(1));
if ~model.massNormalizationSIVerified
    warning('pg_fe_ring_modal_model:UnverifiedSIModalMass', ...
        ['The asset modal shapes are not verified as unit-mass normalized ' ...
        'after SI conversion; absolute acceleration and FRF magnitudes ' ...
        'are exploratory.']);
end
if isfield(a,'sensor_mode_shapes_exact') && logical(a.sensor_mode_shapes_exact(1))
    model.sensorObservationStatus = arrayfun(@(angleDeg) sprintf( ...
        'exact exported SENSOR_%.9g patch',angleDeg), ...
        model.sensorAnglesLocalDeg,'UniformOutput',false);
else
    % Preserve the legacy asset interpretation.  Older MAT files were built
    % by rotating/calibrating SENSOR_0 and must not be reported as exact FE
    % sensor observations.
    model.sensorObservationStatus = arrayfun(@(angleDeg) sprintf( ...
        '%.9g-deg FE circumferential proxy',angleDeg), ...
        model.sensorAnglesLocalDeg,'UniformOutput',false);
    if ~isempty(model.sensorObservationStatus)
        model.sensorObservationStatus{1} = 'exact exported SENSOR_0 patch';
    end
end
if size(model.inputShapeNormal,2) ~= numel(model.frequencyHz) || ...
        size(model.sensorShapeRadial,2) ~= numel(model.frequencyHz)
    error('pg_fe_ring_modal_model:ModeCount','Mode counts do not agree.');
end
if size(model.inputShapeNormal,1) ~= numel(model.contactAngleRad)
    error('pg_fe_ring_modal_model:ContactCount', ...
        'Contact angle and input-shape row counts do not agree.');
end
if size(model.sensorShapeRadial,1) ~= numel(model.sensorAnglesLocalDeg) || ...
        numel(model.sensorAnglesGlobalDeg) ~= numel(model.sensorAnglesLocalDeg)
    error('pg_fe_ring_modal_model:SensorCount', ...
        'Sensor angle and radial-shape row counts do not agree.');
end
if any(~isfinite(model.frequencyHz)) || any(model.frequencyHz<=0) || ...
        any(~isfinite(model.inputShapeNormal),'all') || ...
        any(~isfinite(model.sensorShapeRadial),'all')
    error('pg_fe_ring_modal_model:Finite', ...
        'Modal frequencies and mode-shape coefficients must be finite.');
end
end
