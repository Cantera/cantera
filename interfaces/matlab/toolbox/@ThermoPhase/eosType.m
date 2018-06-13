function e = eosType(tp)
% EOSTYPE  Get the type of the equation of state.
% e = eosType(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     An string identifying the equation of state.
%

e = phase_get(tp.tp_id, 43);
