function z = grid(s, name)
% GRID - the grid in one domain.
%

n = domainIndex(s, name);
d = s.domains(n);
z = gridPoints(d);
