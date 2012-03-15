function v = thermalConductivity(a)
% THERMALCONDUCTIVITY  Thermal conductivity in W/m^2/K.
v = trans_get(a.id, 2);
if v == -1.0
    error(geterr);
elseif v < 0.0
    error('exception raised');
end
