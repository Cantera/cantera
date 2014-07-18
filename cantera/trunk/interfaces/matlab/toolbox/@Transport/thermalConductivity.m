function v = thermalConductivity(a)
% THERMALCONDUCTIVITY  Get the thermal conductivity.
% v = thermalConductivity(a)
% :param a:
%     Instance of class :mat:func:`Transport` (or another
%     object derived from Transport)
%     for which the thermal conductivity is desired.
% :return:
%     Thermal conductivity. Units: W/m-K
%

v = trans_get(a.id, 2);
if v == -1.0
    error(geterr);
elseif v < 0.0
    error('exception raised');
end
