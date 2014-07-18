function mmm = meanMolarMass(tp)
% MEANMOLARMASS - Mean molar mass [kg/kmol].
%
%   The mean molar mass is the mole-fraction-weighted sum of the
%   molar masses of the individual species in the phase.
%

warning('Deprecated in favor of meanMolecularWeight.m')
mmm = meanMolecularWeight(tp);
