function v = electricalConductivity(a)
% ELECTRICALCONDUCTIVITY  Electrical conductivity in S/m.
v = trans_get(a.id, 3);
if v == -1.0
    error(geterr);
elseif v < 0.0
    error('exception raised');
end
