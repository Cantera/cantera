function v = intEnergy_mole(a)
% INTENERGY_MOLE - Molar internal energy [J/kmol].
v = thermo_get(a.tp_id,3);
