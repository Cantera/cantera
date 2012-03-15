function i = domainIndex(d)
% DOMAININDEX - domain index.
%
%     This function returns an integer flag denoting the location
%     of the domain, beginning with 1 at the left.
%
i = domain_methods(d.dom_id, 13) + 1;

