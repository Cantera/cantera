function setState_ST(tp, st)
% SETSTATE_ST  Set the specific entropy and temperature.
% setState_ST(tp,st)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param st:
%     Vector of length 2 containing the desired values for the specific
%     entropy (J/kg-K) and temperature (K).
%

thermo_set(tp.tp_id, 27, st);
