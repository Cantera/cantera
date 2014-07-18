function k = speciesIndex(tp, name)
% SPECIESINDEX  Get the index of a species given the name.
% k = speciesIndex(tp,name)
% The index is an integer assigned to each species in sequence as it
% is read in from the input file.
%
% NOTE: In keeping with the conventions used by Matlab, this method
% returns 1 for the first species, 2 for the second, etc. In
% contrast, the corresponding method speciesIndex in the Cantera C++
% and Python interfaces returns 0 for the first species, 1 for the
% second one, etc. ::
%
%     >> ich4 = speciesIndex(gas, 'CH4');
%     >> iho2 = speciesIndex(gas, 'HO2');
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param name:
%     If name is a single string, the return value will be a integer
%     containing the corresponding index. If it is an cell array of
%     strings, the output will be an array of the same shape
%     containing the indices.
% :return:
%     Scalar or array of integers
%

if iscell(name)
    [m,n] = size(name);
    k = zeros(m, n);
    for i = 1:m
        for j = 1:n
            k(i,j) = phase_get(tp.tp_id, 12, name{i,j});
        end
    end
elseif ischar(name)
    k = phase_get(tp.tp_id, 12, name);
else
    error('name must be either a string or cell array of strings')
end
