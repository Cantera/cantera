function n = molarDensity(p)
% MOLARDENSITY - Molar density [kmol/m^3].
%
n = phase_get(p.tp_id,3);
