function b = isothermalCompressibility(a)
% ISOTHERMALCOMPRESSIBILITY - Isothermal Compressibility [1/Pa]
%
%   b = isothermalCompressibility(a)
%       Return the isothermal compressibility of ThermoPhase a
%       in units of 1/Pa
%
b = thermo_get(a.tp_id,26);