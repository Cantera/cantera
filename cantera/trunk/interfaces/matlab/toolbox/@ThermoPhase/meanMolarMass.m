function wtm = meanMolarMass(tp)
% MEANMOLARMASS - Mean molar mass [kg/kmol].
%
%   The mean molar mass is the mole-fraction-weighted sum of the
%   molar masses of the individual species in the phase.
%

wtm = phase_get(tp.tp_id, 4);
