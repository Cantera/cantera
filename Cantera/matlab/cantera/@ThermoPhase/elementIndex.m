function n = elementIndex(a,name)
% ELEMENTINDEX -  The element index of the element with name
% 'name'. 
% 
%   The index is an integer assigned to each element in sequence as it
%   is read in from the input file.
%
%   NOTE: In keeping with the conventions used by Matlab, this method
%   returns 1 for the first element. In contrast, the corresponding
%   method elementIndex in the Cantera C++ and Python interfaces
%   returns 0 for the first element, 1 for the second one, etc.
%
%      ic = elementIndex(gas, 'C');
%      ih = elementIndex(gas, 'H');
%
  n = phase_get(a.tp_id,13,name);