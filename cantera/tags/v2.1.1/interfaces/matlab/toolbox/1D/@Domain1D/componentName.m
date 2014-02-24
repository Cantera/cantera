function s = componentName(d, n)
% COMPONENTNAME - Name of component n.
%
m = length(n);
s = cell(m);
for i = 1:m
    s{i} = domain_methods(d.dom_id, 40, n(i));
end
