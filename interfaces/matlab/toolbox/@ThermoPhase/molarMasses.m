function mm = molarMasses(tp)
% MOLARMASSES - Array of species molar masses [kg/kmol].
%

warning('Deprecated in favor of molecularWeights.m')
mm = molecularWeights(tp);
