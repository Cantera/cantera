function rho = density(tp)
% DENSITY - Mass density [kg/m^3].
%

rho = phase_get(tp.tp_id, 2);
