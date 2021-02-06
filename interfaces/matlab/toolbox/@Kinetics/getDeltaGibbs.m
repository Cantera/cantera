function dG = getDeltaGibbs(a)
% GETDELTAGIBBS  Get the Gibbs free energy of reaction for each reaction.
% dG = getDeltaGibbs(a)
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics) for which the Gibbs free
%     energies of reaction are desired.
% :return:
%     Returns a vector of the Gibbs free energy of reaction
%     for each reaction. Units: J/kmol
%

dG = kinetics_get(a.id, 17, 1);
