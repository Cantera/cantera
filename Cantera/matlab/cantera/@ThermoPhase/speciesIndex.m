function n = speciesIndex(a,name)
% SPECIESINDEX -  The species index of species with name 'name'. 
%
%   The index is an integer assigned to each species in sequence as it
%   is read in from the input file.
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
n = phase_get(a.tp_id,12,name);