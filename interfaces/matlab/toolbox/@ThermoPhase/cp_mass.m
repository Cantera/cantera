function v = cp_mass(tp)
% CP_MASS - Specific heat at constant pressure [J/kg-K].

v = thermo_get(tp.tp_id, 13);
