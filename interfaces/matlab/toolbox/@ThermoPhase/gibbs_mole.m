function v = gibbs_mole(tp)
% GIBBS_MOLE  Get the mole specific Gibbs function.
% v = gibbs_mole(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of mole specific Gibbs functions of the species.
%     Units: J/kmol
%

v = thermo_get(tp.tp_id, 5);
