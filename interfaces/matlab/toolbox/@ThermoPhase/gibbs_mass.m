function v = gibbs_mass(tp)
% GIBBS_MASS  Get the mass specific Gibbs function.
% v = gibbs_mass(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Mass specific Gibbs function of the mixture. Units: J/kg
%

v = thermo_get(tp.tp_id, 12);
