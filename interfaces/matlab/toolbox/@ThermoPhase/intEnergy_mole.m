function v = intEnergy_mole(tp)
% INTENERGY_MOLE  Get the mole specific internal energy.
% v = intEnergy_mole(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Molar specific internal energy of the mixture. Units: J/kmol
%

v = thermo_get(tp.tp_id, 3);
