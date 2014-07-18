function setThermalConductivity(tr, lam)
% SETTHERMALCONDUCTIVITY  Set the thermal conductivity.
% setThermalConductivity(tr, lam)
% This method can only be used with transport models that
% support directly setting the value of the thermal
% conductivity.
%
% :param tr:
%     Instance of class :mat:func:`Transport` (or another
%     object derived from Transport)
% :param lam:
%     Thermal conductivity in W/(m-K)
%

setParameters(tr, 1, 0, lam);
