function v = satTemperature(a, p)
% SATTEMPERATURE - Saturation temperature for pressure p.
v = thermo_get(a.tp_id,23,p);
