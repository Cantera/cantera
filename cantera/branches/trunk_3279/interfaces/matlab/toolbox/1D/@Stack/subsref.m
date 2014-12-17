function b = subsref(s, index)
% SUBSREF  Redefine subscripted references.
% b = subsref(s,index)
% :param s:
%     Instance of class :mat:func:`Stack`
% :param index:
% :return:
%

switch index.type
    case '()'
        b = s.domains(index.subs{:});
    case '.'
        n = domainIndex(s, index.subs);
        b = s.domains(n);
    otherwise
        error('syntax error');
end
