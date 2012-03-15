function v = cv_mass(a)
% CV_MASS - Specific heat at constant volume [J/kg-K].
v = thermo_get(a.tp_id,14);
