function nsp = nTotalSpecies(a)
% NTOTALSPECIES  Get the total number of species.
% nsp = nTotalSpecies(a)
% The total number of species, summed over all
% participating phases.
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which the number of species is desired.
% :return:
%     Integer total number of species
%

nsp = kinetics_get(a.id, 3, 0);
