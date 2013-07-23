function y = massFraction(r, species)
% MASSFRACTION - Mass fraction of species with name 'species'.
%
k = speciesIndex(r.contents, species) - 1;
y = reactormethods(30, reactor_hndl(r), k);
