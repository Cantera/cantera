function a = setState_satVapor(tp)
% SETSTATE_SATVAPOR    Set the fluid to the saturated vapor state at the
% current temperature.
%

thermo_set(tp.tp_id, 3, 0);
