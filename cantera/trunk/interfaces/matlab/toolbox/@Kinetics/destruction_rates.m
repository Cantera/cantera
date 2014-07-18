function ddot = destruction_rates(a)
% destruction_rates  Chemical destruction rates for all species.
%
%    q = destruction_rates(a)
%
%        Returns a column vector of the destruction rates of all species.
%
%    See also: creation_rates, net_production_rates.
%

ddot = destructionRates(a);
if nargout == 0
    figure
    set(gcf, 'Name', 'Destruction Rates')
    bar(ddot)
    xlabel('Species Number')
    ylabel('Destruction Rate (kmol/m^3/s)')
    title('Species Chemical Destruction Rates')
end
