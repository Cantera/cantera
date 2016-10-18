function clear(tp)
% CLEAR  Delete the kernel object.
% clear(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
%

thermo_set(tp.tp_id, 10, 0);

