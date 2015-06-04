function v = refPressure(tp)
% REFPRESSURE  Get the reference pressure.
% v = refPressure(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Reference pressure in Pa. All standard-state
%     thermodynamic properties are for this pressure.
%

v = thermo_get(tp.tp_id, 15);
