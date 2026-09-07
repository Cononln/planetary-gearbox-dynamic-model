function [force,detail] = pg_coupled_force_loaded_v2(t,p,meta)
%PG_COUPLED_FORCE_LOADED_V2 TE forcing plus physical drive/load torques.

[force,detail.meshExcitation] = pg_coupled_force(t,p,meta);
[operatingForce,detail.operating] = ...
    pg_operating_force_v2(t,p,numel(force));
force = force+operatingForce;
end
