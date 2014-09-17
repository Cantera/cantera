function v = enthalpy_mass(a)
% ENTHALPY_MASS - Specific enthalpy [J/kg].
%
v = thermo_get(a.tp_id,9);
