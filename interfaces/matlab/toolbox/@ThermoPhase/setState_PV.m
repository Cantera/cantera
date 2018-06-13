function setState_PV(tp, pv)
% SETSTATE_PV  Set the pressure and specific volume.
% setState_PV(tp,pv)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param pv:
%     Vector of length 2 containing the desired values for the
%     pressure (Pa) and specific volume (m^3/kg).
%

thermo_set(tp.tp_id, 29, pv);
