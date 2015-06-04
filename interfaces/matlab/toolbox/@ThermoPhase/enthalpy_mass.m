function v = enthalpy_mass(tp)
% ENTHALPY_MASS  Get the mass specific enthalpy.
% v = enthalpy_mass(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of mass specific enthalpies of the species. Units: J/kg
%

v = thermo_get(tp.tp_id, 9);
