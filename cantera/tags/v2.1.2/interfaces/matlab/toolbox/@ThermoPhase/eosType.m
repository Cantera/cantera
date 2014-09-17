function e = eosType(a)
% EOSTYPE - Equation of state type.
%
%    This method returns an integer flag identifying the type of
%    equation of state.
%
e = thermo_get(a.tp_id, 18);
