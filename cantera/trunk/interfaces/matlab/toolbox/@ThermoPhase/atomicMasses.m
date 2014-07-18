function x = atomicMasses(tp)
% ATOMICMASSES  Get the atomic masses of the elements.
% x = atomicMasses(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase).
% :return:
%     Vector of element atomic masses. Units: kg/kmol
%

x = phase_get(tp.tp_id, 30);
