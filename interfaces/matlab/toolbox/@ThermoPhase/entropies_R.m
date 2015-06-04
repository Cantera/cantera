function s = entropies_R(tp)
% ENTROPIES_R  Get the non-dimensional entropy.
% s = entropies_R(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of species non-dimensional entropies.
%

s = thermo_get(tp.tp_id, 36);
