function b = subsref(a, s)
% SUBSREF  Redefine subscripted references for functors.
% b = subsref(a, s)
% :param a:
%     Instance of class :mat:func:`Func`
% :param s:
%     Value at which the function should be evaluated.
% :return:
%     Returns the value of the function evaluated at ``s``
%

switch s.type
    case '()'
        ind = s.subs{:};
        b = zeros(1, length(ind));
        for k = 1:length(ind)
            b(k) = funcmethods(2, a.index, ind(k));
        end
    otherwise
        error('Specify value for x as p(x)')
end
