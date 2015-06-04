function kr = revRateConstants(a)
% REVRATECONSTANTS  Get the reverse reaction rate constants.
% kr = revRateConstants(a)
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
