function v = maxTemp(tp)
% MAXTEMP  Get the maximum temperature of the parameter fits.
% v = maxTemp(tp)
% The parameterizations used to represent the temperature-dependent
% species thermodynamic properties are generally only valid in some
% finite temperature range, which may be different for each species
% in the phase.
%
% See also: :mat:func:`minTemp`
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of maximum temperatures of all species
%

v = thermo_get(tp.tp_id, 17);
