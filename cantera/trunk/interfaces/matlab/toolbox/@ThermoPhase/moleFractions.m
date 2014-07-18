function x = moleFractions(tp)
% MOLEFRACTIONS  Get the mole fractions of all species.
% x = moleFractions(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     Vector of species mole fractions for input phase. If
%     no output argument is specified, a bar plot is produced.
%

x = phase_get(tp.tp_id, 20);
if nargout == 0
    figure
    set(gcf, 'Name', 'Mole Fractions')
    bar(x)
    xlabel('Species Number')
    ylabel('Mole Fraction')
    title('Species Mole Fractions')
end
