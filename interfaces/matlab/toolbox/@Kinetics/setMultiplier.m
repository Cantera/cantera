function setMultiplier(a, irxn, v)
% SETMULTIPLIER  Set the multiplier for the reaction rate of progress.
% setMultiplier(a,irxn,v)
% The multiplier multiplies the reaction rate of progress. It may
% be used to implement sensitivity analysis, or to selectively
% disable reactions.  For reversible reactions, it multiplies both
% the forward and reverse rates. By default, the multiplier value
% is 1.0, but the current value may be checked by calling method
% :mat:func:`multiplier`.
%
% If only two arguments are given, it is assumed that the second is
% the desired multiplication factor for all of the reactions.
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which the multipliers should be set.
% :param irxn:
%     Integer or vector of integers. Reaction number(s) for which
%     the multiplier should be set. Optional.
% :param v:
%     Value by which the reaction rate of progress should be multiplied
%

if nargin == 2
    v = irxn;
    m = nReactions(a);
    irxn = (1:m)';
    n = 1;
elseif nargin == 3
    [m, n] = size(irxn);
else
    error('setMultiplier requires 2 or 3 arguments.')
end

for jm = 1:m
    for jn = 1:n
        kinetics_set(a.id, 1, irxn(jm,jn), v);
    end
end
