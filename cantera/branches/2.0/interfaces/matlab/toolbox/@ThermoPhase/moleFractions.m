function x = moleFractions(a)
% MOLEFRACTIONS - Mole fractions.
%
%     moleFractions(phase)
%
%   returns the array of species mole fractions for phase 'phase'.  If
%   no output argument is specified, a bar plot is produced.
%
x = phase_get(a.tp_id,20);
if nargout == 0
    figure
    set(gcf,'Name','Mole Fractions')
    bar(x)
    xlabel('Species Number')
    ylabel('Mole Fraction')
    title('Species Mole Fractions')
end
