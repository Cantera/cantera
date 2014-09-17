function v = vaporFraction(a)
% VAPORFRACTION - Vapor fraction.
%   If object 'a' represents a liquid/vapor mixture, this method
%   returns the vapor fraction.
v = thermo_get(a.tp_id,22);
