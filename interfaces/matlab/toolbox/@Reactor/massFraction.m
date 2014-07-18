function y = massFraction(r, species)
% MASSFRACTION - Mass fraction of species with name 'species'.
%

if ischar(species)
    k = speciesIndex(r.contents, species) - 1;
else
    k = species - 1;
end

y = reactormethods(30, reactor_hndl(r), k);
