function ddot = destruction_rates(a)
% DESTRUCTION_RATES  Get the chemical destruction rates.
% ddot = destruction_rates(a)
% This function is deprecated in favor of the function
% :mat:func:`destructionRates`
%

warning('This function is deprecated. Use destructionRates instead.')
ddot = destructionRates(a);
if nargout == 0
    figure
    set(gcf, 'Name', 'Destruction Rates')
    bar(ddot)
    xlabel('Species Number')
    ylabel('Destruction Rate (kmol/m^3/s)')
    title('Species Chemical Destruction Rates')
end
