function dS = getDeltaEntropies(a)
% GETDELTAENTROPIES  Get the entropy of reaction for each reaction.
% dS = getDeltaEntropies(a)
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics) for which the entropies of
%     reaction are desired.
% :return:
%     Returns a vector of the entropy of reaction for each
%     reaction. Units: J/kmol-K
%

dS = kinetics_get(a.id, 17, 2);
