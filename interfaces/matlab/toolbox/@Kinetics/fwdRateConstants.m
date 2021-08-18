function kf = fwdRateConstants(a)
% FWDRATECONSTANTS  Get the forward reaction rate constants.
% kf = fwdRateConstants(a)
%
% The computed values include all temperature-dependent, pressure-dependent, and
% third body contributions. Units are a combination of kmol, m^3 and s, that
% depend on the rate expression for the reaction.
%
% Deprecated 2.6:
%
%     Behavior to change after Cantera 2.6; for Cantera 2.6, rate constants of
%     three-body reactions are multiplied with third-body concentrations
%     (no change to legacy behavior). After Cantera 2.6, results will no longer
%     include third-body concentrations and be consistent with conventional
%     definitions (see Eq. 9.75 in Kee, Coltrin and Glarborg, 'Chemically
%     Reacting Flow', Wiley Interscience, 2003).
%     To switch to new behavior, run 'useLegacyRateConstants(0)'.
%
% see also: :mat:func:`revRateConstants`, :mat:func:`equil_Kc`
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which forward rate constants are desired.
% :return:
%     Returns a column vector of the forward rate constants of
%     all of the reactions.
%

kf = kinetics_get(a.id, 15, 0);
