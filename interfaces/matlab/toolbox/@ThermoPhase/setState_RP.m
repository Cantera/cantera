function setState_RP(tp, rp)
% SETSTATE_RP  Set the density and pressure.
% setState_RP(tp, [density,p])
% The density is set first, then the pressure is set by
% changing the temperature holding the density and
% chemical composition fixed.
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param rp:
%     Vector of length 2 containing the desired values for the density (kg/m^3)
%     and pressure (Pa)
%
warning('To be removed after Cantera 3.0. Renamed to setState_DP.')
thermo_set(tp.tp_id, 26, rp);
