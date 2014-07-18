function v = critPressure(tp)
% CRITPRESSURE - Critical pressure [Pa].
%

v = thermo_get(tp.tp_id, 20);
