function Psi = pg_fe_ring_input_shape(theta,model)
%PG_FE_RING_INPUT_SHAPE Periodic interpolation of FE contact participation.
theta = mod(theta(:),2*pi);
a = model.contactAngleRad(:);
v = model.inputShapeNormal;
Psi = interp1([a(end)-2*pi;a;a(1)+2*pi], ...
    [v(end,:);v;v(1,:)],theta,'linear');
if any(~isfinite(Psi),'all')
    error('pg_fe_ring_input_shape:Interpolation','Non-finite interpolation.');
end
end
