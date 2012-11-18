function a = isInlet(d)
% ISINLET - Returns 1 if the domain is an inlet, and 0 otherwise.
%   
t = domainType(d);
if t == 104
  a = 1;
else
  a = 0;
end
