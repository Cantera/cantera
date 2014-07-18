function n = componentIndex(d, name)
% COMPONENTINDEX  Get the index of a component given its name.
% n = componentIndex(d, name)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :param name:
%     String name of the component to look up. If a numeric value
%     is passed, it will be returned.
% :return:
%     Index of the component, or input numeric value.
%

if isa(name, 'double')
    n = name;
else
    n = domain_methods(d.dom_id, 18, name);
end
