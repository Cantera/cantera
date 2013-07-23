function a = setState_Psat(a,px)
% SETSTATE_PSAT    Set the fluid to a saturated state at
% pressure P
%
thermo_set(a.tp_id,24,px);
