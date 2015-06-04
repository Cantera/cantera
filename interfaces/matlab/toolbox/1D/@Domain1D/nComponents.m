function n = nComponents(d)
% NCOMPONENTS  Get the number of components.
% n = nComponents(d)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :return:
%     Number of variables at each grid point
%

n = domain_methods(d.dom_id, 11);
