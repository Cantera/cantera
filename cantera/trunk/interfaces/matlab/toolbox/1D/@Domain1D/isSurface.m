function a = isSurface(d)
% ISSURFACE - Returns 1 if the domain is a surface, and 0 otherwise.
%
t = domainType(d);
if t == 102
    a = 1;
else
    a = 0;
end
