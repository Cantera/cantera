function v = entropy_mole(tp)
% ENTROPY_MOLE  Get the mole specific entropy.
% v = entropy_mole(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of molar specific entropies of the species. Units: J/kg
%

v = thermo_get(tp.tp_id, 4);
