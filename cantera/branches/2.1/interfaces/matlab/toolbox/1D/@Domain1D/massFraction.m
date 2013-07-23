function y = massFraction(d, k)
% MASSFRACTION - Mass fraction of species k.
%
%    This method returns the mass fraction of species k, where
%    k is the integer index of the species in the flow domain
%    to which the boundary domain is attached.
%
if domainIndex(d) == 0
    error('no flow domain attached!')
end

if isInlet(d)
    y = domain_methods(d.dom_id,16,k-1);
else
    error('not yet...');
end
