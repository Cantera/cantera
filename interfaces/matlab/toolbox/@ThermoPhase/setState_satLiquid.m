function setState_satLiquid(tp)
% SETSTATE_SATLIQUID  Set the fluid to the saturated liquid state at the current temperature.
% setState_satLiquid(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
%

thermo_set(tp.tp_id, 2, 0);
