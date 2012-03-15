function x = massFractions(a)
% MASSFRACTIONS - Mass fractions.
%
%     massFractions(phase);
%
%   returns the array of species mass fractions for phase 'phase'.  If
%   no output argument is specified, a bar plot is produced.
%
x = phase_get(a.tp_id,21);
if nargout == 0
    figure
    set(gcf,'Name','Mass Fractions')
    bar(x)
    xlabel('Species Number')
    ylabel('Mass Fraction')
    title('Species Mass Fractions')
end
