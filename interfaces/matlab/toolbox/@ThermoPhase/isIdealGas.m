function v = isIdealGas(tp)
% ISIDEALGAS  Get a flag indicating whether the phase is an ideal gas.
% v = isIdealGas(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     True (1) if the phase is an ideal gas or ideal gas
%     mixture, and false (0) otherwise.
%

if strcmp(eosType(tp), 'ideal-gas')
    v = 1;
else
    v = 0;
end
