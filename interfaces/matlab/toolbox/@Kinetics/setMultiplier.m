function setMultiplier(a, irxn, v)
% SETMULTIPLIER  Set the rate of progress multiplier.
%
%    SETMULTIPLIER(K, IRXN, V) sets the multipler for reaction IRXN
%    to value V.
%
%    see also: MULTIPLIER
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
