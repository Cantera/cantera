function n = domainIndex(d, name)
% DOMAININDEX - Index of the domain with a specified name.
if isa(name,'double')
    n = name;
else
    n = stack_methods(d.stack_id, 109, name);
end
