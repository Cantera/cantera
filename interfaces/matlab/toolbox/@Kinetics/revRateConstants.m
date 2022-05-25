function kr = revRateConstants(a)
% REVRATECONSTANTS  Get the reverse reaction rate constants.
% kr = revRateConstants(a)
%
% The computed values include all temperature-dependent and pressure-dependent
% contributions. By default, third-body concentrations are only considered if
% they are part of the reaction rate definition; for a legacy implementation that
% includes third-body concentrations see :mat:func:`useLegacyRateConstants`.
% Units are a combination of kmol, m^3 and s, that depend on the rate expression
% for the reaction.
%
% See also: :mat:func:`fwdRateConstants`, :mat:func:`equil_KC`
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which reverse rate constants are desired.
% :return:
%     Returns a column vector of the reverse rate constants of
%     all of the reactions.
%

kr = kinetics_get(a.id, 16, 0);
