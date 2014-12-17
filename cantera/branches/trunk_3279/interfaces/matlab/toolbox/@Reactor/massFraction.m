function y = massFraction(r, species)
% MASSFRACTION  Get the mass fraction of a species.
% y = massFraction(r, species)
% :param r:
%     Instance of class :mat:func:`Reactor`
% :param species:
%     String or the one-based integer index of the species
% :return:
%     The mass fraction of the specifed species in the reactor at the
%     end of the last call to :mat:func:`advance` or :mat:func:`step`.
%

if ischar(species)
    k = speciesIndex(r.contents, species) - 1;
else
    k = species - 1;
end

y = reactormethods(30, reactor_hndl(r), k);
