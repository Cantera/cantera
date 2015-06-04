function mmm = meanMolarMass(tp)
% MEANMOLARMASS  Get the mean molecular weight.
% wtm = meanMolarMass(tp)
% Deprecated in favor of :mat:func:`meanMolecularWeight`
%
% See also: :mat:func:`meanMolecularWeight`
%

warning('Deprecated in favor of meanMolecularWeight.m')
mmm = meanMolecularWeight(tp);
