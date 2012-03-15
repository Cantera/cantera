function setCoverages(s,cov)
% SETCOVERAGES - set surface coverages
%
if isa(cov,'double')
    sz = length(cov);
    if sz == nSpecies(s)
        surfmethods(thermo_hndl(s), 3, cov);
    else
        error('wrong size for coverage array');
    end
elseif isa(cov,'char')
    surfmethods(thermo_hndl(s), 5, cov);
end
