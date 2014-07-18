function setState_satLiquid(tp)
% SETSTATE_SATLIQUID    Set the fluid to the saturated liquid state
% at the current temperature.
%

thermo_set(tp.tp_id, 2, 0);
