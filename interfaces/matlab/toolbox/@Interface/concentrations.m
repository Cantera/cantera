function c = concentrations(s)
% CONCENTRATIONS  Get the concentrations of the species on an interface.
% c = concentrations(s)
% :param s:
%     Instance of class :mat:func:`Interface` with surface species
% :return:
%     If no output value is assigned, a bar graph will be plotted.
%     Otherwise, a vector of length ``n_surf_species`` will be
%     returned.
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
