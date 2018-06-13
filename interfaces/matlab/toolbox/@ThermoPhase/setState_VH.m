function setState_VH(tp, vh)
% SETSTATE_VH  Set the specific volume and specific enthalpy.
% setState_VH(tp,vh)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param vh:
%     Vector of length 2 containing the desired values for the specific
%     volume (m^3/kg) and specific enthalpy (J/kg).
%

thermo_set(tp.tp_id, 31, vh);
