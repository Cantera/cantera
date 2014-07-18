function mm = molarMasses(tp)
% MOLARMASSES  Get the molecular weights of all the species.
% x = molarMasses(a)
% Deprecated in favor of :mat:func:`molecularWeights`
%
% See also: :mat:func:`molecularWeights`
%

warning('Deprecated in favor of molecularWeights.m')
mm = molecularWeights(tp);
