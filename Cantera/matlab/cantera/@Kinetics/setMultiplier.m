function setMultiplier(a,i,v)
% SETMULTIPLIER  Set the rate of progress multiplier.
%
%    SETMULTIPLIER(K, IRXN, V) sets the multipler for reaction IRXN
%    to value V.
%
%    see also: MULTIPLIER
%
kinetics_set(a.id,1,i,v);

