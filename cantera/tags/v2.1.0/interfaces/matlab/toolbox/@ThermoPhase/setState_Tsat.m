function a = setState_Tsat(a,tx)
% SETSTATE_TSAT    Set the fluid to a saturated state at
% temperature t
%
thermo_set(a.tp_id,25,tx);
