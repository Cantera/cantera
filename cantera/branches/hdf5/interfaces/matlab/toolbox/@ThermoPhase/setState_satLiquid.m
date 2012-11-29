function a = setState_satLiquid(a)
% SETSTATE_SATLIQUID    Set the fluid to the saturated liquid state
% at the current temperature.
%
thermo_set(a.tp_id,2,0);
