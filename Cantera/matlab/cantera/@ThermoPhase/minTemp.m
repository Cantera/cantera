function v = critTemperature(a)
% CRITTEMPERATURE - Critical temperature [K].
%
%    The critical temperature is the temperature at the critical point
%
  v = thermo_get(p.tp_id,19);
