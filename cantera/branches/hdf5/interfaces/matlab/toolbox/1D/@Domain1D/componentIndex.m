function n = componentIndex(d, name)
% COMPONENTINDEX -
%
if isa(name,'double')
    n = name;
else
    n = domain_methods(d.dom_id, 18, name);
end
