function cdot = creationRates(a)
% CREATIONRATES  Chemical creation rates (kmol/m^3/s).
%
%    cdot = creationRates(K)
%
%        Returns a column vector of the creation rates of all
%        species. If the output is not assigned to a variable, a
%        bar graph is produced.
%
%    See also: destructionRates, netProdRates.
%
cdot = kinetics_get(a.id,21,0);
if nargout == 0
    figure
    set(gcf,'Name','Creation Rates')
    bar(cdot)
    xlabel('Species Number')
    ylabel('Creation Rate (kmol/m^3-s)')
    title('Species Chemical Creation Rates')
end
