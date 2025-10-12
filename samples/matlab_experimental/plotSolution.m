function plotSolution(domain, component)
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

    z = domain.grid;
    x = domain.values(component);
    plot(z, x);
    xlabel('z (m)');
    ylabel(component);
end
