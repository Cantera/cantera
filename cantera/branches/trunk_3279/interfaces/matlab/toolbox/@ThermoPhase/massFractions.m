function y = massFractions(tp)
% MASSFRACTIONS  Get the mass fractions of all species.
% y = massFractions(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of species mass fractions for input phase. If
%     no output argument is specified, a bar plot is produced.
%

y = phase_get(tp.tp_id, 21);
if nargout == 0
    figure
    set(gcf, 'Name', 'Mass Fractions')
    bar(y)
    xlabel('Species Number')
    ylabel('Mass Fraction')
    title('Species Mass Fractions')
end
