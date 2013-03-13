function v = satPressure(a, T)
% SATPRESSURE - Saturation pressure for temperature T.
v = thermo_get(a.tp_id,24,T);
