function v = useLegacyRateConstants(legacy)
% USELEGACYRATECONSTANTS  Set definition used for rate constant calculation
% useLegacyRateConstants(0)
%
% If set to 0 (false - default value), rate constants of three-body reactions are
% consistent with conventional definitions (for example Eq. 9.75 in Kee, Coltrin
% and Glarborg, 'Chemically Reacting Flow', Wiley Interscience, 2003). If set to
% 1 (true), output for rate constants of three-body reactions is multiplied by
% third-body concentrations, consistent with Cantera's behavior prior to version 3.0.
%
% See also: :mat:func:`fwdRateConstants`, :mat:func:`revRateConstants`.
%

ctmethods(0, 8, legacy);
