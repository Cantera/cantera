function v = gibbs_mole(tp)
% GIBBS_MOLE - Molar Gibbs function [J/kmol].

v = thermo_get(tp.tp_id, 5);
