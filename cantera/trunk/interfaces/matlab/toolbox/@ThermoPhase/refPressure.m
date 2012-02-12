function v = refPressure(p)
% REFPRESSURE - Reference pressure [Pa].
%
%    All standard-state thermodynamic properties are for this pressure.
%
v = thermo_get(p.tp_id,15);
