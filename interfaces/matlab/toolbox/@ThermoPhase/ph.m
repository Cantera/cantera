function i = ph(tp)
% PH  Get the ph field of the phase.
% i = ph(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%

warning('method ph is deprecated.');
i = tp.ph;
