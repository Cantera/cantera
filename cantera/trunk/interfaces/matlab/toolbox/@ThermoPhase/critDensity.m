function v = critDensity(a)
% CRITDENSITY - Critical density [kg/m3].
%
v = thermo_get(a.tp_id,21);
