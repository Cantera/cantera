function c = concentrations(s)
% CONCENTRATIONS - Surface concentrations
%
c = surfmethods(thermo_hndl(s), 103);
if nargout == 0
    figure
    set(gcf, 'Name', 'Concentrations')
    bar(c);
    colormap(summer);
    nm = speciesNames(s);
    set(gca, 'XTickLabel', nm);
    xlabel('Species Name');
    ylabel('Concentration [kmol/m2]');
    title('Surface Species Concentrations');
end
