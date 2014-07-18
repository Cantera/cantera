function v = vaporFraction(tp)
% VAPORFRACTION - Vapor fraction.
%   If object 'a' represents a liquid/vapor mixture, this method
%   returns the vapor fraction.

v = thermo_get(tp.tp_id, 22);
