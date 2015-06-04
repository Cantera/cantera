function v = satTemperature(tp, p)
% SATTEMPERATURE  Get the saturation temperature for a given pressure.
% v = satTemperature(tp, p)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :param p:
%     Pressure. Units: Pa
% :return:
%     Saturation temperature for pressure p. Units: K
%

v = thermo_get(tp.tp_id, 23, p);
