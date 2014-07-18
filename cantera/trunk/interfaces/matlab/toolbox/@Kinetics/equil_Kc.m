function kc = equil_Kc(a)
% EQUIL_KC  Get the equilibrium constants for all reactions
% kc = equil_Kc(a)
%
% See also: :mat:func:`fwdRateConstants`, :mat:func:`revRateConstants`
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which equilibrium constants are desired.
% :return:
%     Returns a column vector of the equilibrium constants
%     for all reactions. The vector has an entry for every
%     reaction, whether reversible or not, but non-zero values
%     occur only for the reversible reactions. If the output is
%     not assigned to a variable, a bar graph is produced instead.
%

kc = kinetics_get(a.id, 14, 0);
if nargout == 0
    figure
    set(gcf, 'Name', 'Equilibrium Constants')
    bar(log10(kc))
    xlabel('Reaction Number')
    ylabel('log_{10} Kc [kmol, m, s]')
    title('Equilibrium Constants Kc')
end
