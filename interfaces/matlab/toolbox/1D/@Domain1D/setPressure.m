function setPressure(d, p)
% SETPRESSURE  Set the pressure.
% d = setPressure(d, p)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :param p:
%     Pressure to be set. Units: Pa
%

domain_methods(d.dom_id, 63, p);
