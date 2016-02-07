function setCoverages(s, cov, norm)
% SETCOVERAGES  Set surface coverages of the species on an interface.
% setCoverages(s, cov, norm)
% :param s:
%      Instance of class :mat:func:`Interface`
% :param cov:
%      Coverage of the species. ``cov`` can be either a vector of
%      length ``n_surf_species``, or a string in the format
%      ``'Species:Coverage, Species:Coverage'``
% :param norm:
%      Optional flag that denotes whether or not to normalize the species
%      coverages. ``norm`` is either of the two strings ``'nonorm'``` or
%      ``'norm'``. If left unset, the default is `norm`. This only works if
%      ``s`` is a vector, not a string. Since unnormalized coverages can lead to
%      unphysical results, ``'nonorm'`` should be used only in rare cases, such
%      as computing partial derivatives with respect to a species coverage.
%

if nargin == 3 && strcmp(norm,'nonorm')
    norm_flag = 0;
else
    norm_flag  = 1;
end

if isa(cov, 'double')
    sz = length(cov);
    if sz == nSpecies(s)
        surfmethods(thermo_hndl(s), 3, cov, norm_flag);
    else
        error('wrong size for coverage array');
    end
elseif isa(cov,'char')
    surfmethods(thermo_hndl(s), 5, cov);
end
