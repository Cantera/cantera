function wdot = netProdRates(a)
% NETPRODRATES  Get the net chemical production rates for all species.
% wdot = netProdRates(a)
%
% See also: :mat:func:`creationRates`, :mat:func:`destructionRates`
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which net production rates are desired.
% :return:
%     Returns a column vector of the net production (creation -
%     destruction) rates of all species. If the output is not
%     assigned to a variable, a bar plot is produced.
%

wdot = kinetics_get(a.id, 23, 0);
if nargout == 0
    figure
    set(gcf, 'Name', 'Production Rates')
    bar(wdot)
    xlabel('Species Number')
    ylabel('Net Production Rate (kmol/m^3/s)')
    title('Species Net Chemical Production Rates')
end
