function plotSolution(s, domain, component)
    %% Utility for plotting a specific solution component
    %
    %     >> plotSolution(s, domain, component)
    %
    % :s:
    %    Instance of class Matlab class `Sim1D`.
    % :domain:
    %    Name of domain from which the component should be retrieved.
    % :component:
    %    Name of the component to be plotted

    n = s.stackIndex(domain);
    d = s.domains{n};
    z = d.gridPoints;
    x = s.getSolution(domain, component);
    plot(z, x);
    xlabel('z (m)');
    ylabel(component);
end
