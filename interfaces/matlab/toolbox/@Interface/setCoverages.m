function setCoverages(s, cov)
% SETCOVERAGES  Set surface coverages of the species on an interface.
% setCoverages(s, cov)
% :param s:
%      Instance of class :mat:func:`Interface`
% :param cov:
%      Coverage of the species. ``cov`` can be either a vector of
%      length ``n_surf_species``, or a string in the format
%      ``'Species:Coverage, Species:Coverage'``
%

if isa(cov, 'double')
    sz = length(cov);
    if sz == nSpecies(s)
        surfmethods(thermo_hndl(s), 3, cov);
    else
        error('wrong size for coverage array');
    end
elseif isa(cov,'char')
    surfmethods(thermo_hndl(s), 5, cov);
end
