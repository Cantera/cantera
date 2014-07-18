function v = gibbs_mass(tp)
% GIBBS_MASS - Specific Gibbs function [J/kg].

v = thermo_get(tp.tp_id, 12);
