function setState_TV(tp, tv)
% SETSTATE_TV  Set the temperature and specific volume.
% setState_TV(tp,tv)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param tv:
%     Vector of length 2 containing the desired values for the
%     temperature (K) and specific volume (m^3/kg).
%

thermo_set(tp.tp_id, 28, tv);
