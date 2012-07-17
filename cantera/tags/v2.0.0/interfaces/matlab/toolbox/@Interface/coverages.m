function c = coverages(s)
% COVERAGES - Surface coverages
%
c = surfmethods(thermo_hndl(s), 101);
if nargout == 0
    figure
    set(gcf,'Name','Coverages')
    bar(c);
    colormap(summer);
    nm = speciesNames(s);
    legend(nm);
    xlabel('Species Number');
    ylabel('Coverage');
    title('Surface Species Coverages');
end
