function v = cv_mass(tp)
% CV_MASS - Specific heat at constant volume [J/kg-K].

v = thermo_get(tp.tp_id, 14);
