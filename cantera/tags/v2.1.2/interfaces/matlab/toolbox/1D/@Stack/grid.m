function z = grid(s, d)
% GRID - the grid in one domain.
%
n = domainIndex(s,d);
d = s.domains(n);
z = gridPoints(d);
