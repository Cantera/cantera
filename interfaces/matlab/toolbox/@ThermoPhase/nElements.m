function n = nElements(tp)
% NELEMENTS  Get the number of elements.
% n = nElements(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Number of elements in the phase.
%

n = phase_get(tp.tp_id, 10);
