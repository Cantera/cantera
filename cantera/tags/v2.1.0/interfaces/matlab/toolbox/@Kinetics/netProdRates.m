function wdot = netProdRates(a)
% NETPRODRATES  Net chemical production rates for all species.
%
%    wdot = netProdRates(a)
%
%        Returns a column vector of the net production (creation -
%        destruction) rates of all species. If the output is not
%        assigned to a variable, a bar plot is produced.
%
%    See also: creationRates, destructionRates
%
wdot = kinetics_get(a.id,23,0);
if nargout == 0
    figure
    set(gcf,'Name','Production Rates')
    bar(wdot)
    xlabel('Species Number')
    ylabel('Net Production Rate (kmol/m^3/s)')
    title('Species Net Chemical Production Rates')
end
