function v = domainType(d)
% DOMAINTYPE  Get the type of domain.
% v = domainType(d)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :return:
%     This function returns a string describing the domain type.
%

v = domain_methods(d.dom_id, 12);
