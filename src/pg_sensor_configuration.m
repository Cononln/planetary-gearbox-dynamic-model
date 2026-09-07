function sensor = pg_sensor_configuration(p)
%PG_SENSOR_CONFIGURATION Validate fixed sensor positions and sensitive axes.

if isfield(p.sensor,'locationAngles')
    location = p.sensor.locationAngles(:).';
else
    location = p.sensor.angles(:).';
end
if isfield(p.sensor,'sensitiveAngles')
    sensitive = p.sensor.sensitiveAngles(:).';
else
    sensitive = location;
end
if numel(location) ~= numel(sensitive)
    error('pg_sensor_configuration:SizeMismatch', ...
        'Sensor locationAngles and sensitiveAngles must have equal length.');
end
if isempty(location) || any(~isfinite(location)) || any(~isfinite(sensitive))
    error('pg_sensor_configuration:InvalidAngle', ...
        'Sensor angles must be nonempty finite values in radians.');
end

sensor.locationAngles = location;
sensor.sensitiveAngles = sensitive;
sensor.nSensor = numel(location);
sensor.locationAnglesDeg = rad2deg(location);
sensor.sensitiveAnglesDeg = rad2deg(sensitive);
end
