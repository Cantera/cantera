function k = elementIndex(tp, name)
% ELEMENTINDEX  Get the index of an element given its name.
% k = elementIndex(tp,name)
% The index is an integer assigned to each element in sequence as it
% is read in from the input file.
%
% If ``name`` is a single string, the return value will be a integer
% containing the corresponding index. If it is an cell array of
% strings, the output will be an array of the same shape
% containing the indices.
%
% NOTE: In keeping with the conventions used by Matlab, this method
% returns 1 for the first element. In contrast, the corresponding
% method elementIndex in the Cantera C++ and Python interfaces
% returns 0 for the first element, 1 for the second one, etc. ::
%
%     >> ic = elementIndex(gas, 'C');
%     >> ih = elementIndex(gas, 'H');
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :param name:
%     String or cell array of strings of elements to look up
% :return:
%     Integer or vector of integers of element indices
%

if iscell(name)
    [m, n] = size(name);
    k = zeros(m, n);
    for i = 1:m
        for j = 1:n
            k(i,j) = phase_get(tp.tp_id, 13, name{i,j});
        end
    end
else
    k = phase_get(tp.tp_id, 13, name);
end
