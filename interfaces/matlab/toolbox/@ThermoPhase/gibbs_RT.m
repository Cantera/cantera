function g_RT = gibbs_RT(tp)
% GIBBS_RT - Species non-dimensional Gibbs free energies.
%
%        This method returns an array containing the pure species
%        standard-state Gibbs free energies.
%

g_RT = enthalpies_RT(tp) - entropies_R(tp);
