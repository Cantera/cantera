function v = critDensity(tp)
% CRITDENSITY - Critical density [kg/m3].
%

v = thermo_get(tp.tp_id, 21);
