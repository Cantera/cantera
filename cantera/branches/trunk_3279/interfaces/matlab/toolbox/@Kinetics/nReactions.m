function n = nReactions(a)
% NREACTIONS  Get the number of reactions.
% n = nReactions(a)
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which the number of reactions is desired.
% :return:
%     Integer number of reactions
%

n = kinetics_get(a.id, 1, 0);
