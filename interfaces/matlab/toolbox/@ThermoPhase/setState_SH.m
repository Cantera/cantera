function setState_SH(tp, sh)
% SETSTATE_SH  Set the specific entropy and specific enthalpy.
% setState_SH(tp,sh)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param sh:
%     Vector of length 2 containing the desired values for the specific
%     entropy (J/kg/K) and specific enthalpy (J/kg).
%

thermo_set(tp.tp_id, 33, sh);
