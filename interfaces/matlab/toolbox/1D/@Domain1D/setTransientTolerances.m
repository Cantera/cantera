function d = setSteadyTolerances(d, component, rtol, atol)
% SETSTEADYTOLERANCES -
%

if strcmp(component, 'default')
    nc = nComponents(d);
    for ii = 1:nc
        domain_methods(d.dom_id, 56, ii, rtol, atol);
    end
elseif iscell(component)
    nc = length(component);
    for ii = 1:nc
        n = componentIndex(d, component{ii});
        domain_methods(d.dom_id, 56, n, rtol, atol);
    end
else
    n = componentIndex(d, component);
    domain_methods(d.dom_id, 56, n, rtol, atol);
end
