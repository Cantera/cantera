function nu = stoich_net(a,species,rxns)
% stoich_net  Net stoichiometric coefficients.
%
%    nu = stoich_net(a)
%
%        Returns a sparse matrix of all net (product - reactant)
%        stoichiometric coefficients. The matrix element nu(k,i) is the
%        net stoichiometric coefficient of species k in reaction i.
%
%    nu = stoich_net(a, species, rxns)
%
%        Returns a sparse matrix the same size as above, but
%        containing only entries for the specified species and
%        reactions. For example, stoich_net(a,3,[1 3 5 7]) returns a
%        sparse matrix containing only the coefficients for species 3
%        in reactions 1, 3, 5, and 7.
%
%    Note that the net stoichiometric coefficients may be negative,
%    unlike the reactant or product stoichiometric coefficients.
%
%    See also: stoich_r, stoich_p.
%
if nargin == 1
    nu = stoich_p(a) - stoich_r(a);
elseif nargin == 3
    nu = stoich_p(a,species,rxns) - stoich_r(a,species,rxns);
else
    error(['syntax error. Type ''help stoich_net'' for more' ...
        ' information.'])
end
