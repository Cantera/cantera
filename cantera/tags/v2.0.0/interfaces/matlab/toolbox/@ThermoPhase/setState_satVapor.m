function a = setState_satVapor(a)
% SETSTATE_SATVAPOR    Set the fluid to the saturated vapor state at the
% current temperature.
%
thermo_set(a.tp_id,3,0);
