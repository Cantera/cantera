function v = entropy_mass(tp)
% ENTROPY_MASS  Get the mass specific entropy.
% v = entropy_mass(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of mass specific entropies of the species. Units: J/kg-K
%

v = thermo_get(tp.tp_id, 11);
