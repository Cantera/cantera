function v = satTemperature(tp, p)
% SATTEMPERATURE - Saturation temperature for pressure p.

v = thermo_get(tp.tp_id, 23, p);
