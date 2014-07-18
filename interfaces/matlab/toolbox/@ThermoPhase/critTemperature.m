function v = critTemperature(tp)
% CRITTEMPERATURE - Critical temperature [K].
%

v = thermo_get(tp.tp_id, 19);
