function v = cv_mole(a)
% CV_MOLE - Molar heat capacity at constant volume [J/kmol-K].
v = thermo_get(a.tp_id,7);
