function v = useLegacyRateConstants(legacy)
% USELEGACYRATECONSTANTS  Set definition used for rate constant calculation
% useLegacyRateConstants(0)
%
% If set to 0 (false), rate constants of three-body reactions are consistent with
% conventional definitions. If set to 1 (true), output for rate constants of
% three-body reactions is multipied by third-body concentrations (legacy behavior).
% For the pre-compiled Cantera 2.6 distribution, the default value is set to 1 (true),
% which implies no change compared to previous behavior. For user-compiled Cantera,
% the default behavior can be changed by the SCons flag 'legacy_rate_constants'.
%
% Deprecated 2.6:
%
%     Behavior to change after Cantera 2.6; for Cantera 2.6, rate constants of
%     three-body reactions are multiplied with third-body concentrations
%     (no change to legacy behavior). After Cantera 2.6, results will no longer
%     include third-body concentrations and be consistent with conventional
%     definitions (see Eq. 9.75 in Kee, Coltrin and Glarborg, 'Chemically
%     Reacting Flow', Wiley Interscience, 2003).
%
% See also: :mat:func:`fwdRateConstants`, :mat:func:`revRateConstants`.
%

ctmethods(0, 8, legacy);
