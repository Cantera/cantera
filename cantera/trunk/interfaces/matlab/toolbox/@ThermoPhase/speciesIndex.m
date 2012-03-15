function k = speciesIndex(a,name)
% SPECIESINDEX -  The species index of species with name 'name'.
%
%   The index is an integer assigned to each species in sequence as it
%   is read in from the input file.
%
%   If name is a single string, the return value will be a integer
%   containing the corresponding index. If it is an cell array of
%   strings, the output will be an array of the same shape
%   containing the indices.
%
%   NOTE: In keeping with the conventions used by Matlab, this method
%   returns 1 for the first species, 2 for the second, etc. In
%   contrast, the corresponding method speciesIndex in the Cantera C++
%   and Python interfaces returns 0 for the first species, 1 for the
%   second one, etc.
%
%      ich4 = speciesIndex(gas, 'CH4');
%      iho2 = speciesIndex(gas, 'HO2');
%

if iscell(name)
    [m, n] = size(name);
    k = zeros(m,n);
    for i = 1:m
        for j = 1:n
            k(i,j) = phase_get(a.tp_id,12,name{i,j});
        end
    end
else
    k = phase_get(a.tp_id,12,name);
end
