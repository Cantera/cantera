function ddot = destructionRates(a)
% destructionRates  Chemical destruction rates (kmol/m^3/s).
%
%    cdot = destructionRates(a)
%
%        Returns a column vector of the destruction rates of all
%        species. If the output is not assigned to a variable, a
%        bar graph is produced.
%
%    See also: creationRates, netProdRates.
%
ddot = kinetics_get(a.id,22,0);
if nargout == 0
    figure
    set(gcf,'Name','Destruction Rates')
    bar(ddot)
    xlabel('Species Number')
    ylabel('Destruction Rate (kmol/m^3/s)')
    title('Species Chemical Destruction Rates')
end
