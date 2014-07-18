function d = setProfile(d, n, p)
% SETPROFILE -
%
if d.stack == 0
    error('Install domain in stack before calling setProfile.');
end

setProfile(d.stack,domainIndex(d), n, p);
