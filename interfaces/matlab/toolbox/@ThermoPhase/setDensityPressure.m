function setDensityPressure(tp, density, p)
% SETDENSITYPRESSURE  Set the density and pressure.
% setDensityPressure(tp,density,p)
% The density is set first, then the pressure is set by
% changing the temperature holding the density and
% chemical composition fixed.
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param density:
%     Density. Units: kg/m^3
% :param p:
%     Pressure. Units: Pa
%

if p <= 0.0:
    error('The pressure must be positive.')
end

if density <= 0.0:
    error('The pressure must be positive.')
end

thermo_set(tp.tp_id, 26, density, p);
