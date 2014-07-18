function n = multiplier(a,irxn)
% MULTIPLIER  Get the multiplier for reaction rate of progress.
% n = multiplier(a,irxn)
% The multiplier multiplies the reaction rate of progress. It may
% be used to implement sensitivity analysis, or to selectively
% disable reactions.  For reversible reactions, it multiplies both
% the forward and reverse rates. By default, the multiplier value
% is 1.0, but it may be set to any other value by calling method
% :mat:func:`setMultiplier`.
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which the multipliers are desired.
% :param irxn:
%     Integer reaction number for which the multiplier is desired.
% :return:
%     Multiplier of the rate of progress of reaction number ``irxn``
%

n = kinetics_get(a.id, 2, irxn);
