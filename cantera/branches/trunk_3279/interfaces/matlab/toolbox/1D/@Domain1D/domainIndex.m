function i = domainIndex(d)
% DOMAININDEX  Get the domain index.
% i = domainIndex(d)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :return:
%     This function returns an integer flag denoting the location
%     of the domain, beginning with 1 at the left.
%

i = domain_methods(d.dom_id, 13) + 1;

