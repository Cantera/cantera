function y = massFraction(d, k)
% MASSFRACTION  Get the mass fraction of a species given its integer index.
% y = massFraction(d, k)
% This method returns the mass fraction of species ``k``, where
% k is the integer index of the species in the flow domain
% to which the boundary domain is attached.
%
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :param k:
%     Integer species index
% :return:
%     Mass fraction of species
%

if domainIndex(d) == 0
    error('No flow domain attached!')
end

if isInlet(d)
    y = domain_methods(d.dom_id, 16, k-1);
else
    error('Input domain must be an inlet');
end
