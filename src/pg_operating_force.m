function [force,detail] = pg_operating_force(t,p,nSecondOrder)
%PG_OPERATING_FORCE Balanced sun drive and carrier resisting torques.

force = zeros(nSecondOrder,1);
inputTorque = 0;
outputTorque = 0;
rippleFraction = 0;
rippleFrequencyHz = p.kin.f_s;
if isfield(p.operating,'inputTorque'); inputTorque = p.operating.inputTorque; end
if isfield(p.operating,'outputTorque'); outputTorque = p.operating.outputTorque; end
if isfield(p.operating,'torqueRippleFraction')
    rippleFraction = p.operating.torqueRippleFraction;
end
if isfield(p.operating,'torqueRippleFrequencyHz')
    rippleFrequencyHz = p.operating.torqueRippleFrequencyHz;
end

ripple = rippleFraction*cos(2*pi*rippleFrequencyHz*t);
instantaneousInputTorque = inputTorque*(1+ripple);
force(p.map.sun(3)) = instantaneousInputTorque;
force(p.map.carrier(3)) = -outputTorque;

detail.inputTorque = instantaneousInputTorque;
detail.nominalInputTorque = inputTorque;
detail.outputTorque = outputTorque;
detail.ringReactionTarget = outputTorque-instantaneousInputTorque;
detail.rippleFraction = rippleFraction;
end
