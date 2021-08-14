function v = useLegacyRateConstants(legacy)
% USELEGACYRATECONSTANTS  Set definition used for rate constant calculation
% useLegacyRateConstants(1)
%
% If set to 1 (true), rate constants include third-body concentrations for
% ThreeBodyReaction objects.
%
% See also: :mat:func:`fwdRateConstants`, :mat:func:`revRateConstants`.
%

ctmethods(0, 8, legacy);
