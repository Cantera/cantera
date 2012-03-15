function setValveCoeff(f, k)
% SETVALVECOEFF - set valve coefficient
%
if f.type ~= 3
    error('Valve coefficient can only be set for valves')
end
ok = flowdevicemethods(4, f.index, k);
if ok < 0
    error(geterr);
end
