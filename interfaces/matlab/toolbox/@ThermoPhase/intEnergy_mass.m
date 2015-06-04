function v = intEnergy_mass(tp)
% INTENERGY_MASS  Get the mass specific internal energy.
% v = intEnergy_mass(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of mass specific internal energies of the species.
%     Units: J/kg
%

v = thermo_get(tp.tp_id, 10);
