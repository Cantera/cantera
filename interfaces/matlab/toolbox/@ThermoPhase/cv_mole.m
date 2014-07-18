function v = cv_mole(tp)
% CV_MOLE - Molar heat capacity at constant volume [J/kmol-K].

v = thermo_get(tp.tp_id, 7);
