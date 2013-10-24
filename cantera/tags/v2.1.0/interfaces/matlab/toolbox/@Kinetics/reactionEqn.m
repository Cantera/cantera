function e = reactionEqn(a, irxn)
% reactionEqn  Reaction equation of reaction irxn.
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
    e = kinetics_get(a.id, 31, irxn); % rxnstring(a.id, irxn);
else
    e = cell(m,n);
    for i = 1:m
        for j = 1:n
            e{i,j} = kinetics_get(a.id, 31, irxn(i,j)); % rxnstring(a.id, irxn(i,j));
        end
    end
end
