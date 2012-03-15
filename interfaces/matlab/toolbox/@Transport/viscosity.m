function v = viscosity(a)
v = trans_get(a.id, 1);
if v == -1.0
    error(geterr);
elseif v < 0.0
    error('exception raised');
end
