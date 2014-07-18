function v = critPressure(tp)
% CRITPRESSURE  Get the critical pressure.
% v = critPressure(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Critical pressure. Units: Pa
%

v = thermo_get(tp.tp_id, 20);
