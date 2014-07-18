function v = enthalpy_mole(a)
% ENTHALPY_MOLE  Get the mole specific enthalpy.
% v = enthalpy_mole(a)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of molar specific enthalpies of the species. Units: J/kmol
%

v = thermo_get(a.tp_id, 2);
