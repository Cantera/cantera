function setState_DP(tp, dp)
% SETSTATE_DP  Set the density and pressure.
% setState_DP(tp, [density,p])
% The density is set first, then the pressure is set by
% changing the temperature holding the density and
% chemical composition fixed.
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param dp:
%     Vector of length 2 containing the desired values for the density (kg/m^3)
%     and pressure (Pa)
%

thermo_set(tp.tp_id, 26, dp);
