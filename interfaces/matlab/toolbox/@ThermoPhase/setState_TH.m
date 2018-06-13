function setState_TH(tp, th)
% SETSTATE_TH  Set the temperature and specific enthalpy.
% setState_TH(tp,th)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param th:
%     Vector of length 2 containing the desired values for the
%     temperature (K) and specific enthalpy (J/kg).
%

thermo_set(tp.tp_id, 32, th);
