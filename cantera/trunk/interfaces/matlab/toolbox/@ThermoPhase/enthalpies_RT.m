function v = enthalpies_RT(tp)
% ENTHALPIES_RT - Pure species non-dimensional enthalpies.
%
%      h_rt = enthalpies_RT(phase)
%
%    sets array h_rt to the array of standard-state species enthalpies
%    for phase 'phase' divided by RT, where R is the universal gas
%    constant and T is the temperature.  For gaseous species, these
%    values are ideal gas enthalpies.
%
%

v = thermo_get(tp.tp_id, 32);
