function v = cp_mass(a)
% CP_MASS - Specific heat at constant pressure [J/kg-K].
v = thermo_get(a.tp_id,13);
