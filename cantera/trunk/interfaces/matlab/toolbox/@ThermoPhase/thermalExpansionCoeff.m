function a = thermalExpansionCoeff(p)
% THERMALEXPANSIONCOEFF - Thermal Expansion Coefficient [1/K]
%
%   a = thermalExpansionCoeff(p)
%       Return the thermal expansion coefficient of ThermoPhase p
%       in units of 1/K
%
a = thermo_get(p.tp_id,27);
