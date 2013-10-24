function rho = density(p)
% DENSITY - Mass density [kg/m^3].
%
rho = phase_get(p.tp_id,2);
