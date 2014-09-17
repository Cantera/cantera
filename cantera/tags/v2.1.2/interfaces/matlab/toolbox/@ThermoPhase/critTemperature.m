function v = critTemperature(a)
% CRITTEMPERATURE - Critical temperature [K].
%
v = thermo_get(a.tp_id,19);
