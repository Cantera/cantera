function plotSolution(s, domain, component)
    % Plot a specified solution component. ::
    %
    %     >> plotSolution(s, domain, component)
    %
    % :param s:
    %    Instance of class :mat:class:`Sim1D`.
    % :param domain:
    %    Name of domain from which the component should be
    %    retrieved.
    % :param component:
    %    Name of the component to be plotted.

    n = s.stackIndex(domain);
    d = s.domains{n};
    z = d.gridPoints;
    x = s.getSolution(domain, component);
    plot(z, x);
    xlabel('z (m)');
    ylabel(component);
end