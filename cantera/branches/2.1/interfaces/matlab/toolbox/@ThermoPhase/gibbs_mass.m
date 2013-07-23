function v = gibbs_mass(a)
% GIBBS_MASS - Specific Gibbs function [J/kg].
v = thermo_get(a.tp_id,12);
