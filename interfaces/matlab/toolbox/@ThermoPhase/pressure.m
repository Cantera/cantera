function v = pressure(tp)
% PRESSURE  Get the pressure.
% v = pressure(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Pressure. Units: Pa
%

v = thermo_get(tp.tp_id, 8);
