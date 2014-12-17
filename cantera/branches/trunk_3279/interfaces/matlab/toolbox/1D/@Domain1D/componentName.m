function s = componentName(d, n)
% COMPONENTNAME  Get the name of a component given its index.
% s = componentName(d, n)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :param n:
%     Integer or vector of integers of components' names
%     to get.
% :return:
%     Cell array of component names.
%

m = length(n);
s = cell(m);
for i = 1:m
    s{i} = domain_methods(d.dom_id, 40, n(i));
end
