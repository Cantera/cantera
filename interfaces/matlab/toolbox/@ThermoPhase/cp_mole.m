function v = cp_mole(tp)
% CP_MOLE - Molar heat capacity at constant pressure [J/kmol-K].

v = thermo_get(tp.tp_id, 6);
