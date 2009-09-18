function yn = isReversible(a, i)
% ISREVERSIBLE - Reversible reaction flag.
%
%    A reversible reaction is one that runs in both the forward
%    direction (reactants -> products) and in the reverse direction
%    (products -> reactants). The reverse rate for reversible
%    reactions is computed from thermochemistry, so that the
%    reaction satisfies detailed balance, and the net rate of
%    progress is zero in states of chemical equilibrium. 
%
%       ISREVERSIBLE(K, IRXN) returns 1 if reaction number IRXN is 
%       reversible, and 0 if it is irreversible.
%
yn = kinetics_get(a.id,4,i);
