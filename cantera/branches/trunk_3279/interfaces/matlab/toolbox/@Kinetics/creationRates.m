function cdot = creationRates(a)
% CREATIONRATES  Get the chemical creation rates.
% cdot = creationRates(a)
%
% See also: :mat:func:`destructionRates`, :mat:func:`netProdRates`
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which creation rates are desired.
% :return:
%     Returns a column vector of the creation rates of all
%     species. If the output is not assigned to a variable, a
%     bar graph is produced. Units: kmol/m**3-s
%

cdot = kinetics_get(a.id, 21, 0);
if nargout == 0
    figure
    set(gcf, 'Name', 'Creation Rates')
    bar(cdot)
    xlabel('Species Number')
    ylabel('Creation Rate (kmol/m^3-s)')
    title('Species Chemical Creation Rates')
end
