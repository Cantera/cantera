function b = subsref(a,s)
% SUBSREF 
switch s.type
case '()'
    ind = s.subs{:};
    for k = 1:length(ind)
        b(k) = funcmethods(2,a.index,ind(k))
    end
otherwise
   error('Specify value for x as p(x)')
end
