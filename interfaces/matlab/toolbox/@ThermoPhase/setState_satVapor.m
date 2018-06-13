function setState_satVapor(tp)
% SETSTATE_SATVAPOR  Set the fluid to the saturated vapor state at the current temperature.
% setState_satVapor(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
%

set(tp, 'T', temperature(tp), 'Vapor', 1.0)
