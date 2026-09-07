function Psi = pg_ring_mode_shape(theta,model)
%PG_RING_MODE_SHAPE Evaluate radial analytical modes at circumferential angle.
% Rows correspond to theta samples; columns correspond to modes 7 through 12.

theta = theta(:);
Psi = zeros(numel(theta),numel(model.frequencyHz));
for r = 1:numel(model.frequencyHz)
    argument = model.waveNumber(r)*theta;
    if strcmp(model.family{r},'cos')
        Psi(:,r) = cos(argument);
    else
        Psi(:,r) = sin(argument);
    end
end
end
