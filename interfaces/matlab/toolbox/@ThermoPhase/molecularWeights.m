function mw = molecularWeights(tp)
% MOLECULARWEIGHTS  Get the molecular weights of the species.
% x = molecularWeights(tp)
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of species molecular weights. Units: kg/kmol
%

mw = phase_get(tp.tp_id, 22);
