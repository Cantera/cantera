function setState_UP(tp, up)
% SETSTATE_UP  Set the specific internal energy and pressure.
% setState_UP(tp,up)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param up:
%     Vector of length 2 containing the desired values for the specific
%     internal energy (J/kg) and pressure (Pa).
%

thermo_set(tp.tp_id, 30, up);
