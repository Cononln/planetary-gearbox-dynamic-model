function [force,detail] = pg_operating_force_v2(t,p,nSecondOrder)
%PG_OPERATING_FORCE_V2 Apply balanced drive and load torques.
% Retained as a compatibility wrapper for the loaded-v2 workflow.
[force,detail] = pg_operating_force(t,p,nSecondOrder);
end
