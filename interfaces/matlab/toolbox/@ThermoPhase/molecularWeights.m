function mw = molecularWeights(tp)
% MOLECULARWEIGHTS - Array of species molar masses [kg/kmol].
%
%   This method is deprecated - use molarMasses instead.
%
%   See also: molarMasses
%

mw = phase_get(tp.tp_id, 22);
