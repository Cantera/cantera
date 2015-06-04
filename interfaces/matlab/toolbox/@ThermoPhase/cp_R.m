function v = cp_R(tp)
% CP_R  Get the non-dimensional specific heats at constant pressure.
% v = cp_R(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of specific heats of the species at
%     constant pressure, non-dimensional basis
%

v = thermo_get(tp.tp_id, 38);
