function v = cv_mass(tp)
% CV_MASS  Get the mass-basis specific heats at constant volume.
% v = cv_mass(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of specific heats of the species at
%     constant volume. Units: J/kg-K
%

v = thermo_get(tp.tp_id, 14);
