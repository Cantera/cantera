function c = coverages(s)
% COVERAGES  Get the surface coverages of the species on an interface.
% c = coverages(s)
% :param s:
%     Instance of class :mat:func:`Interface` with surface species
% :return:
%     If no output value is assigned, a bar graph will be plotted.
%     Otherwise, a vector of length ``n_surf_species`` will be
%     returned.
%

c = surfmethods(thermo_hndl(s), 101);
if nargout == 0
    figure
    set(gcf, 'Name', 'Coverages')
    bar(c);
    colormap(summer);
    nm = speciesNames(s);
    set(gca,'XTickLabel', nm);
    xlabel('Species Name');
    ylabel('Coverage');
    title('Surface Species Coverages');
end
