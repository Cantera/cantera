function zz = gridPoints(d, n)
% GRIDPOINTS  Get grid points from a domain.
% zz = gridPoints(d, n)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :param n:
%     Optional, vector of grid points to be retrieved.
% :return:
%     Vector of grid points. Length of ``n`` or :mat:func:`nPoints`.
%

if nargin == 1
    zz = zeros(1, nPoints(d));
    for i = 1:nPoints(d)
        zz(i) = domain_methods(d.dom_id, 19, i);
    end
else
    m = length(n);
    zz = zeros(1, m);
    for i = 1:m
        zz(i) = domain_methods(d.dom_id, 19, n(i));
    end
end
