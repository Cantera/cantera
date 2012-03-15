function v = gibbs_mole(a)
% GIBBS_MOLE - Molar Gibbs function [J/kmol].
v = thermo_get(a.tp_id,5);
