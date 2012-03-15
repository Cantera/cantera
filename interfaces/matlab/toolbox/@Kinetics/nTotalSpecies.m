function nsp = nTotalSpecies(a)
% NTOTALSPECIES - The total number of species, summed over all
% participating phases.
%
nsp = kinetics_get(a.id, 3, 0);
