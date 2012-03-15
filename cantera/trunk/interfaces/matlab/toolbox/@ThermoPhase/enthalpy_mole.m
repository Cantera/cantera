function v = enthalpy_mole(a)
% ENTHALPY_MOLE - Molar enthalpy [J/kmol].
%
v = thermo_get(a.tp_id,2);
