function v = critPressure(a)
% CRITPRESSURE - Critical pressure [Pa].
%
v = thermo_get(a.tp_id,20);
