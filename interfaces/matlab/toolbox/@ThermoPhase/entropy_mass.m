function v = entropy_mass(tp)
% ENTROPY_MASS  Get the mass specific entropy.
% v = entropy_mass(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Mass specific entropy of the mixture. Units: J/kg-K
%

v = thermo_get(tp.tp_id, 11);
