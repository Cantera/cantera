function setState_SP(tp, sp)
% SETSTATE_SP    Set the specific entropy [J/kg/K] and pressure [Pa].
%
%    setState_SP(a, sp) sets the specific entropy and pressure
%    of object a, holding its composition fixed. Argument 'sp' must
%    be a vector of length 2 containing the desired values for the specific
%    entropy (J/kg/K) and pressure (Pa).
%

if sp(1) <= 0.0
    error('The specific entropy must be positive.');
end
if sp(2) <= 0.0
    error('The pressure must be positive.');
end
thermo_set(tp.tp_id, 23, sp);
