function [M,C,K,meta] = pg_assemble_pengyue_system(t,p,fault)
%PG_ASSEMBLE_PENGYUE_SYSTEM Carrier-fixed 18-DOF reproduction assembly.
% Peng Yue defines all perturbation coordinates in the OXY frame fixed to
% and rotating with the carrier, then neglects gyroscopic/centrifugal terms
% at the reported low speed. Consequently the mesh/pin direction vectors are
% fixed in that frame while the TVMS and transmission error remain periodic.

[stiffness,tvmsDetail] = pg_tvms(t,p,fault);
[M,C,K,meta] = pg_assemble_baseline(0,p,stiffness);
meta.tvms = tvmsDetail;
meta.fault = fault;
meta.coordinateFrame = 'carrier-fixed rotating frame; G and Komega neglected';
end
