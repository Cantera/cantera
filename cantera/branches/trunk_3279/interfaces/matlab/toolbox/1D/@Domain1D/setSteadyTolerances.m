function setSteadyTolerances(d, component, rtol, atol)
% SETSTEADYTOLERANCES  Set the steady-state tolerances.
% d = setSteadyTolerances(d, component, rtol, atol)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :param component:
%     String or cell array of strings of component values
%     whose tolerances should be set. If ``'default'`` is
%     specified, the tolerance of all components will be set.
% :param rtol:
%     Relative tolerance
% :param atol:
%     Absolute tolerance
%

if strcmp(component, 'default')
    nc = nComponents(d);
    for ii = 1:nc
        domain_methods(d.dom_id, 55, ii, rtol, atol);
    end
elseif iscell(component)
    nc = length(component);
    for ii = 1:nc
        n = componentIndex(d, component{ii});
        domain_methods(d.dom_id, 55, n, rtol, atol);
    end
else
    n = componentIndex(d, component);
    domain_methods(d.dom_id, 55, n, rtol, atol);
end
