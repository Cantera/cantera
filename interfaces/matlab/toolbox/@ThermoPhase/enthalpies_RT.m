function v = enthalpies_RT(tp)
% ENTHALPIES_RT  Get the non-dimensional enthalpies.
% v = enthalpies_RT(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of standard-state species enthalpies
%     divided by RT, where R is the universal gas
%     constant and T is the temperature. For gaseous species, these
%     values are ideal gas enthalpies.
%

v = thermo_get(tp.tp_id, 32);
