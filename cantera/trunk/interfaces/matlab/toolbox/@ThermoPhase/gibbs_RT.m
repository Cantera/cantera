function g_RT = gibbs_RT(p)
% GIBBS_RT - Species non-dimensional Gibbs free energies.
%
%        This method returns an array containing the pure species
%        standard-state Gibbs free energies.
%
g_RT = enthalpies_RT(p) - entropies_R(p);
