function v = entropies_R(p)
% ENTROPIES_R - Species non-dimensional entropies.
%
%        This method returns an array containing the pure species
%        standard-state entropies.
%
v = thermo_get(p.tp_id,36);
