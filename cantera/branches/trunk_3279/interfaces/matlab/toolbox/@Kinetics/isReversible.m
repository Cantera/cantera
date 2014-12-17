function yn = isReversible(a, i)
% ISREVERSIBLE  Get an array of flags indicating reversibility of a reaction.
% yn = isReversible(a, i)
% A reversible reaction is one that runs in both the forward
% direction (reactants -> products) and in the reverse direction
% (products -> reactants). The reverse rate for reversible
% reactions can computed from thermochemistry, so that the
% reaction satisfies detailed balance, and the net rate of
% progress is zero in states of chemical equilibrium. The reverse
% rate can also be specified directly by a rate expression. An
% irreversible reaction is one whose reverse reaction rate is
% zero.
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which the reversible flags are desired.
% :param i:
%     Integer reaction number
% :return:
%     1 if reaction number ``i`` is
%     reversible, and 0 if it is irreversible.
%

yn = kinetics_get(a.id, 4, i);
