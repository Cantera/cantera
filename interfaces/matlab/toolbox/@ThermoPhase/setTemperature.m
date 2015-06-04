function setTemperature(tp,t)
% SETTEMPERATURE  Set the temperature.
% setTemperature(tp,t)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param t:
%     Temperature. Units: K
%

if (t <= 0)
    error('The temperature must be positive.');
end
phase_set(tp.tp_id, 1, t);
