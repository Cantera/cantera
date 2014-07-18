function v = enthalpy_mass(tp)
% ENTHALPY_MASS - Specific enthalpy [J/kg].
%

v = thermo_get(tp.tp_id, 9);
