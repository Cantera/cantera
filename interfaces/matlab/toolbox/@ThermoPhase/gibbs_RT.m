function g_RT = gibbs_RT(tp)
% GIBBS_RT  Get the non-dimensional Gibbs function.
% g_RT = gibbs_RT(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of non-dimensional Gibbs functions of the species.
%

g_RT = enthalpies_RT(tp) - entropies_R(tp);
