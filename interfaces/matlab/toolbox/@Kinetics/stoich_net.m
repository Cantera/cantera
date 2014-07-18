function nu = stoich_net(a, species, rxns)
% STOICH_NET  Get the net stoichiometric coefficients.
% nu = stoich_net(a,species,rxns)
%
% See also: :mat:func:`stoich_r`, :mat:func:`stoich_p`
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which the net stoichiometric coefficients are desired.
% :param species:
%     Species indices for which net stoichiometric coefficients
%     should be retrieved. Optional argument; if specified, ``rxns``
%     must be specified as well.
% :param rxns:
%     Reaction indices for which net stoichiometric coefficients
%     should be retrieved. Optional argument; if specified, ``species``
%     must be specified as well.
% :return:
%     Returns a sparse matrix of all net stoichiometric
%     coefficients. The matrix element ``nu(k,i)`` is the
%     stoichiometric coefficient of species k as a net in
%     reaction i. If ``species`` and ``rxns`` are specified, the matrix
%     will contain only entries for the specified species and
%     reactions. For example, ``stoich_p(a,3,[1 3 5 7])`` returns a
%     sparse matrix containing only the coefficients for species 3
%     in reactions 1, 3, 5, and 7.
%

if nargin == 1
    nu = stoich_p(a) - stoich_r(a);
elseif nargin == 3
    nu = stoich_p(a, species, rxns) - stoich_r(a, species, rxns);
else
    error(['stoich_net requires 1 or 3 arguments.'])
end
