function v = intEnergy_mass(tp)
% INTENERGY_MASS  Get the mass specific internal energy.
% v = intEnergy_mass(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Mass specific internal energy of the mixture. Units: J/kg
%

v = thermo_get(tp.tp_id, 10);
