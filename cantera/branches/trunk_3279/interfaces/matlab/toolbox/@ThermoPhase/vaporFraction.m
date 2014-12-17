function v = vaporFraction(tp)
% VAPORFRACTION  Get the vapor fraction.
% v = vaporFraction(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :return:
%     Vapor fraction.
%

v = thermo_get(tp.tp_id, 22);
