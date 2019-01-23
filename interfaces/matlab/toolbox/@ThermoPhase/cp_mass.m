function v = cp_mass(tp)
% CP_MASS  Get the mass-basis specific heat at constant pressure.
% v = cp_mass(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Mass basis specific heat of the mixture at
%     constant pressure. Units: J/kg-K
%

v = thermo_get(tp.tp_id, 13);
