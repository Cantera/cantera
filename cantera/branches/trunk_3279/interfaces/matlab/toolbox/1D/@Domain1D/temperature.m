function t = temperature(d)
% TEMPERATURE  Get the boundary temperature.
% t = temperature(d)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :return:
%     Temperature. Units: K
%

t = domain_methods(d.dom_id, 15);
