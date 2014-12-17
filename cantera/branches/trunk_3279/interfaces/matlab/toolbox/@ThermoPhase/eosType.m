function e = eosType(tp)
% EOSTYPE  Get the type of the equation of state.
% e = eosType(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     An integer flag identifying the type of equation of state.
%     See the definitions in include/cantera/thermo/mix_defs.h
%

e = thermo_get(tp.tp_id, 18);
