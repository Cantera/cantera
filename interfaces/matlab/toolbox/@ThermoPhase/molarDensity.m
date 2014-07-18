function n = molarDensity(tp)
% MOLARDENSITY - Molar density [kmol/m^3].
%

n = phase_get(tp.tp_id, 3);
