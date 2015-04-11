function clear(d)
% CLEAR  Delete the C++ Domain1D object.
% clear(d)
% :param d:
%     Instance of class :mat:func:`Domain1D` (or another
%     object that derives from Domain1D)
%

domain_methods(d.dom_id, 10);
