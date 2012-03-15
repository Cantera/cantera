function v = intEnergy_mass(a)
% INTENERGY_MASS - Specific internal energy [J/kg].
v = thermo_get(a.tp_id,10);
