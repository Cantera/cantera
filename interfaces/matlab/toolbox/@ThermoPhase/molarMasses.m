function x = molarMasses(tp)
% MOLARMASSES - Array of species molar masses [kg/kmol].
%

x = phase_get(tp.tp_id, 22);
