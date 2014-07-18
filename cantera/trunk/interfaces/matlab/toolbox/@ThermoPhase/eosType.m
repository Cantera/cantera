function e = eosType(tp)
% EOSTYPE - Equation of state type.
%
%    This method returns an integer flag identifying the type of
%    equation of state.
%

e = thermo_get(tp.tp_id, 18);
