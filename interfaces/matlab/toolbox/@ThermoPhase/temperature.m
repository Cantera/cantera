function t = temperature(tp)
% TEMPERATURE  Get the temperature.
% t = temperature(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :return:
%     Temperature. Units: K
%

t = phase_get(tp.tp_id, 1);
