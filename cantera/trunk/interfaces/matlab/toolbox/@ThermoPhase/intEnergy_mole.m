function v = intEnergy_mole(tp)
% INTENERGY_MOLE - Molar internal energy [J/kmol].
v = thermo_get(tp.tp_id, 3);
