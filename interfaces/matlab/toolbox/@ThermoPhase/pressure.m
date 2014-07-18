function v = pressure(tp)
% PRESSURE - Pressure [Pa].
%

v = thermo_get(tp.tp_id, 8);
