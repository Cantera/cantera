function plotSolution(s, domain, component)
% PLOTSOLUTION - plot a specified solution component
%
%     plotSolution(s, 'flow', 'T') plots component 'T' in domain 'flow'
%
n = domainIndex(s,domain);
d = s.domains(n);
z = gridPoints(d);
x = solution(s, domain, component);
plot(z, x);
xlabel('z (m)');
ylabel(component);
