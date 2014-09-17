function b = subsref(s,index)
% SUBSREF -
switch index.type
    case '()'
        b = s.domains(index.subs{:});
    case '.'
        n = domainIndex(s, index.subs);
        b = s.domains(n);
    otherwise
        error('syntax error');
end
