function z = grid(s, name)
% GRID  Get the grid in one domain.
% z = grid(s, name)
% :param s:
%     Instance of class :mat:func:`Stack`
% :param name:
%     Name of the domain for which the grid
%     should be retrieved.
% :return:
%     The grid in domain name
%

n = domainIndex(s, name);
d = s.domains(n);
z = gridPoints(d);
