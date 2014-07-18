function v = intEnergy_mass(tp)
% INTENERGY_MASS - Specific internal energy [J/kg].

v = thermo_get(tp.tp_id, 10);
