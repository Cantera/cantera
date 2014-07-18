function s = entropies_R(tp)
% ENTROPIES_R - Species non-dimensional entropies.
%
%        This method returns an array containing the pure species
%        standard-state entropies.
%

s = thermo_get(tp.tp_id, 36);
