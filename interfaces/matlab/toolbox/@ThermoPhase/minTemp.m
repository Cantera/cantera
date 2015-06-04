function v = minTemp(tp)
% MINTEMP  Get the minimum temperature of the parameter fits.
% v = minTemp(tp)
% The parameterizations used to represent the temperature-dependent
% species thermodynamic properties are generally only valid in some
% finite temperature range, which may be different for each species
% in the phase.
%
% See also: :mat:func:`maxTemp`
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of minimum temperatures of all species
%

v = thermo_get(tp.tp_id, 16);
