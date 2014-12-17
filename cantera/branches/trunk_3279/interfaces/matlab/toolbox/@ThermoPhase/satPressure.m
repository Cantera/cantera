function v = satPressure(tp, T)
% SATPRESSURE  Get the saturation pressure for a given temperature.
% v = satPressure(tp, T)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :param T:
%     Temperature Units: K
% :return:
%     Saturation pressure for temperature T. Units: Pa
%

v = thermo_get(tp.tp_id, 24, T);
