function d = setBounds(d, component, lower, upper)
% SETBOUNDS -
%
n = componentIndex(d,component);
domain_methods(d.dom_id, 51, n, lower, upper);
