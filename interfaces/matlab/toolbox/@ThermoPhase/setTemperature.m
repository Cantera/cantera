function setTemperature(tp,t)
% SETTEMPERATURE  Set the temperature [K].
%

if (t <= 0)
    error('The temperature must be positive.');
end
phase_set(tp.tp_id, 1, t);
