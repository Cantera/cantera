function ddot = destructionRates(a)
% DESTRUCTIONRATES  Get the chemical destruction rates.
% ddot = destructionRates(a)
%
% See also: :mat:func:`creationRates`, :mat:func:`netProdRates`
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which destruction rates are desired.
% :return:
%     Returns a column vector of the destruction rates of all
%     species. If the output is not assigned to a variable, a
%     bar graph is produced. Units: kmol/m**3-s
%

ddot = kinetics_get(a.id, 22, 0);
if nargout == 0
    figure
    set(gcf, 'Name', 'Destruction Rates')
    bar(ddot)
    xlabel('Species Number')
    ylabel('Destruction Rate (kmol/m^3/s)')
    title('Species Chemical Destruction Rates')
end
