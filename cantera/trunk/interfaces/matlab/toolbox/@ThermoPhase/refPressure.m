function v = refPressure(tp)
% REFPRESSURE - Reference pressure [Pa].
%
%    All standard-state thermodynamic properties are for this pressure.
%

v = thermo_get(tp.tp_id, 15);
