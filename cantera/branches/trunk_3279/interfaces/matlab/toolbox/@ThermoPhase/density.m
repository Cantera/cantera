function rho = density(tp)
% DENSITY  Get the density.
% rho = density(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Mass density. Units: kg/m**3
%

rho = phase_get(tp.tp_id, 2);
