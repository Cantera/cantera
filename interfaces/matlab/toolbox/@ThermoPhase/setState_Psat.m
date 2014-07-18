function setState_Psat(tp, px)
% SETSTATE_PSAT    Set the fluid to a saturated state at
% pressure P
%
thermo_set(tp.tp_id, 24, px);
