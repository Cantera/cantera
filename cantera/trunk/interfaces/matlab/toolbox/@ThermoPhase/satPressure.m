function v = satPressure(tp, T)
% SATPRESSURE - Saturation pressure for temperature T.

v = thermo_get(tp.tp_id, 24, T);
