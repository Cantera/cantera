function n = multiplier(a,irxn)
% MULTIPLIER  Multiplier for reaction rate of progress.
%
%    The multiplier multiplies the reaction rate of progress. It may
%    be used to implement sensitivity analysis, or to selectively
%    disable reactions.  For reversible reactions, it multiplies both
%    the forward and reverse rates. By default, the multiplier value
%    is 1.0, but it may be set to any other value by calling method
%    setMultiplier.
%
%       MULTIPLIER(K, IRXN)   Multiplier for reaction number IRXN
%
n = kinetics_get(a.id,2,irxn);
