function v = cp_mole(a)
% CP_MOLE - Molar heat capacity at constant pressure [J/kmol-K].
v = thermo_get(a.tp_id,6);
