function a = setTemperature(a,t)
% SETTEMPERATURE  Set the temperature [K].
%
if (t <= 0)
    error('the temperature must be positive');
end
phase_set(a.tp_id,1,t);
