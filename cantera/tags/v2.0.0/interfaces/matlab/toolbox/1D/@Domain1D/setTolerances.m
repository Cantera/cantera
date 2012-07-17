function d = setTolerances(d, component, rtol, atol, typ)
% SETTOLERANCES -
%
ityp = 0;
if nargin == 5
    switch typ
        case 'ts'
            ityp = -1;
        case 'time'
            ityp = -1;
        case 'ss'
            ityp = 1;
        case 'steady'
            ityp = 1;
    end
end

if strcmp(component,'default')
    nc = nComponents(d);
    for ii = 1:nc
        domain_methods(d.dom_id, 52, ii, rtol, atol, ityp);
    end
    return
end

if iscell(component)
    nc = length(component);
    for ii = 1:nc
        n = componentIndex(d, component{ii});
        domain_methods(d.dom_id, 52, n, rtol, atol, ityp);
    end
else
    n = componentIndex(d, component);
    domain_methods(d.dom_id, 52, n, rtol, atol, ityp);
end
