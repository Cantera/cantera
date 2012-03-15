function d = setProfile(d, n, p)
% SETPROFILE -
%
if d.stack == 0
    error('install domain in stack before calling setProfile.');
end

setProfile(d.stack,domainIndex(d),n,p);
