function a = setState_Tsat(tp, tx)
% SETSTATE_TSAT    Set the fluid to a saturated state at
% temperature t
%

thermo_set(tp.tp_id, 25, tx);
