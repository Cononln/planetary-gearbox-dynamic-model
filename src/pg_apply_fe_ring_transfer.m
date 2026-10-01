function path = pg_apply_fe_ring_transfer(time,meshForceRingPlanet, ...
    carrierAngleRad,p,model,options)
%PG_APPLY_FE_RING_TRANSFER One-way PR force -> FE modes -> sensor response.
% The default observation uses the exported contact and sensor mode shapes.
% Optional kinematic observation operators are retained for sensitivity
% comparisons; they are not calibrated ANSYS contact-to-sensor FRFs.

arguments
    time (:,1) double
    meshForceRingPlanet (:,:) double
    carrierAngleRad (:,1) double
    p struct
    model struct
    options.SensorObservationMethod (1,1) string {mustBeMember( ...
        options.SensorObservationMethod, ...
        ["exported_mode_shapes","rotated_sensor0_kinematic", ...
        "kinematic_hann_modal"])} = ...
        "exported_mode_shapes"
end
time = time(:);
carrierAngleRad = carrierAngleRad(:);
n = numel(time);
if size(meshForceRingPlanet,1) ~= n || ...
        size(meshForceRingPlanet,2) ~= p.model.nPlanet
    error('pg_apply_fe_ring_transfer:ForceSize','PR force size is invalid.');
end
dt = median(diff(time));
if numel(carrierAngleRad) ~= n || any(abs(diff(time)-dt)>1e-8*dt)
    error('pg_apply_fe_ring_transfer:Time','Time/angle vectors are invalid.');
end
fs = 1/dt;
nMode = numel(model.frequencyHz);
nSensor = size(model.sensorShapeRadial,1);
fSigned = (0:n-1)'*fs/n;
fSigned(fSigned>fs/2) = fSigned(fSigned>fs/2)-fs;
w = 2*pi*fSigned;
Y = complex(zeros(n,nSensor));
modalForceRms = zeros(nSensor,nMode);
pathWeights = ones(n,nSensor,p.model.nPlanet);

if options.SensorObservationMethod=="kinematic_hann_modal"
    % Full 18-DOF PR dynamic forces remain the excitation.  Kinematics only
    % supplies a parameter-free periodic transmission-path weight between
    % each moving PR contact and each fixed sensor:
    % h(delta)=0.5*(1+cos(delta)), delta wrapped on the ring.
    % The weighted force then passes through the same co-located SENSOR_0 FE
    % modal kernel for all three sensor angles.  No sensor-specific gain is
    % fitted, so differences come from the actual planet/sensor geometry.
    sensorGlobalRad = deg2rad(model.sensorAnglesGlobalDeg);
    inputAtSensor0 = pg_fe_ring_input_shape(sensorGlobalRad(1),model);
    for is = 1:nSensor
        qModal = zeros(n,nMode);
        for ip = 1:p.model.nPlanet
            thetaGlobal = pg_fe_ring_contact_angle(carrierAngleRad,p,model,ip);
            weight = 0.5*(1+cos(thetaGlobal-sensorGlobalRad(is)));
            pathWeights(:,is,ip) = weight;
            qModal = qModal+(meshForceRingPlanet(:,ip).*weight).* ...
                inputAtSensor0;
        end
        modalForceRms(is,:) = sqrt(mean(qModal.^2,1));
        for im = 1:nMode
            wn = model.omega(im);
            zeta = model.zeta(im);
            Hacc = -(w.^2)./(wn^2-w.^2+1i*2*zeta*wn*w);
            modalAcceleration = Hacc.*fft(qModal(:,im));
            Y(:,is) = Y(:,is)+ ...
                modalAcceleration*model.sensorShapeRadial(1,im);
        end
    end
    methodText = ['18-DOF PR dynamic mesh forces weighted by the exact ' ...
        'planet/sensor angular separation using h=0.5(1+cos(delta)), then ' ...
        'filtered through the co-located mass-normalized SENSOR_0 ANSYS ' ...
        'modal kernel; periodic steady-state modal frequency response'];
elseif options.SensorObservationMethod=="rotated_sensor0_kinematic"
    % SENSOR_0 is the only exact exported observation.  For sensors mounted
    % at known circumferential offsets, rotate the source position into the
    % SENSOR_0 coordinate system and reuse the same physical transfer kernel.
    % This is a circumferential-symmetry assumption, not an empirical gain.
    sensorOffsetRad = deg2rad(model.sensorAnglesGlobalDeg- ...
        model.sensorAnglesGlobalDeg(1));
    for is = 1:nSensor
        qModal = zeros(n,nMode);
        for ip = 1:p.model.nPlanet
            thetaGlobal = pg_fe_ring_contact_angle(carrierAngleRad,p,model,ip);
            thetaEquivalentAtSensor0 = thetaGlobal-sensorOffsetRad(is);
            qModal = qModal+meshForceRingPlanet(:,ip).* ...
                pg_fe_ring_input_shape(thetaEquivalentAtSensor0,model);
        end
        modalForceRms(is,:) = sqrt(mean(qModal.^2,1));
        for im = 1:nMode
            wn = model.omega(im);
            zeta = model.zeta(im);
            Hacc = -(w.^2)./(wn^2-w.^2+1i*2*zeta*wn*w);
            modalAcceleration = Hacc.*fft(qModal(:,im));
            Y(:,is) = Y(:,is)+ ...
                modalAcceleration*model.sensorShapeRadial(1,im);
        end
    end
    methodText = ['one-way PR mesh-force projection through mass-normalized ' ...
        'ANSYS modes; exact SENSOR_0 kernel rotated by kinematic sensor ' ...
        'offset; periodic steady-state modal frequency response'];
else
    qModal = zeros(n,nMode);
    for ip = 1:p.model.nPlanet
        theta = pg_fe_ring_contact_angle(carrierAngleRad,p,model,ip);
        qModal = qModal+meshForceRingPlanet(:,ip).* ...
            pg_fe_ring_input_shape(theta,model);
    end
    modalForceRms(:,:) = repmat(sqrt(mean(qModal.^2,1)),nSensor,1);
    for im = 1:nMode
        wn = model.omega(im);
        zeta = model.zeta(im);
        Hacc = -(w.^2)./(wn^2-w.^2+1i*2*zeta*wn*w);
        modalAcceleration = Hacc.*fft(qModal(:,im));
        Y = Y+modalAcceleration*model.sensorShapeRadial(:,im).';
    end
    methodText = ['one-way PR mesh-force projection through mass-normalized ' ...
        'ANSYS modes; exported/proxy sensor mode shapes; periodic ' ...
        'steady-state modal frequency response'];
end

path.time = time;
path.fs = fs;
path.acceleration = real(ifft(Y,[],1));
path.sensorAnglesLocalDeg = model.sensorAnglesLocalDeg;
path.sensorAnglesGlobalDeg = model.sensorAnglesGlobalDeg;
if options.SensorObservationMethod=="kinematic_hann_modal"
    path.sensorObservationStatus = repmat({ ...
        'kinematic Hann path + exact co-located SENSOR_0 FE kernel'}, ...
        1,nSensor);
elseif options.SensorObservationMethod=="rotated_sensor0_kinematic"
    path.sensorObservationStatus = cell(1,nSensor);
    for is = 1:nSensor
        offsetDeg = mod(model.sensorAnglesGlobalDeg(is)- ...
            model.sensorAnglesGlobalDeg(1),360);
        if is==1
            path.sensorObservationStatus{is} = ...
                'exact exported SENSOR_0 patch';
        else
            path.sensorObservationStatus{is} = sprintf( ...
                'SENSOR_0 transfer kernel rotated by %.9g deg',offsetDeg);
        end
    end
else
    path.sensorObservationStatus = model.sensorObservationStatus;
end
path.sensorObservationMethod = options.SensorObservationMethod;
path.modalGeneralizedForceRms = modalForceRms;
path.pathWeights = pathWeights;
path.method = methodText;
end
