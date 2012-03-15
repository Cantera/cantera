function i = domainType(d)
% DOMAINTYPE - Type of domain.
%
%     This function returns an integer flag denoting the domain
%     type.
i = domain_methods(d.dom_id, 12);

