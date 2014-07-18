function plotSolution(s, domain, component)
% PLOTSOLUTION  Plot a specified solution component.
% plotSolution(s, domain, component)
% :param s:
%     Instance of class :mat:func:`Stack`
% :param domain:
%     Name of domain from which the component should be
%     retrieved
% :param component:
%     Name of the component to be plotted
%

n = domainIndex(s, domain);
d = s.domains(n);
z = gridPoints(d);
x = solution(s, domain, component);
plot(z, x);
xlabel('z (m)');
ylabel(component);
