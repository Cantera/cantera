function x = molarMasses(a)
% MOLARMASSES - Array of species molar masses [kg/kmol].
%
x = phase_get(a.tp_id,22);
