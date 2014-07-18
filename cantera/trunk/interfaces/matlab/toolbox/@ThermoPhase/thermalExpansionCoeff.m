function a = thermalExpansionCoeff(tp)
% THERMALEXPANSIONCOEFF - Thermal Expansion Coefficient [1/K]
%
%   a = thermalExpansionCoeff(p)
%       Return the thermal expansion coefficient of ThermoPhase p
%       in units of 1/K
%

a = thermo_get(tp.tp_id, 27);
