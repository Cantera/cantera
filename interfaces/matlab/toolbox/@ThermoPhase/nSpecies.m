function n = nSpecies(tp)
% NSPECIES  Get the number of species.
% n = nSpecies(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Number of species in the phase.
%

n = phase_get(tp.tp_id, 11);
