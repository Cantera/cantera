function e = reactionEqn(a, irxn)
% REACTIONEQN  Get the reaction equation of a reaction.
% e = reactionEqn(a, irxn)
% If only the first argument
% is given, the reaction equations of all of the reactions are
% returned in a cell array. Otherwise, ``irxn`` must be an integer
% or vector of integers.
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which the reaction equations are desired.
% :param irxn:
%     Optional. Integer or vector of integer reaction numbers.
% :return:
%     String or cell array of strings of the reaction equations.
%

if nargin == 1
    m = nReactions(a);
    n = 1;
    irxn = (1:m)';
elseif nargin == 2
    if isa(irxn,'double')
        [m, n] = size(irxn);
    else
        error('reaction number(s) must be numeric');
    end
end

if m == 1 && n == 1
    e = kinetics_get(a.id, 31, irxn);
else
    e = cell(m,n);
    for i = 1:m
        for j = 1:n
            e{i, j} = kinetics_get(a.id, 31, irxn(i,j));
        end
    end
end
