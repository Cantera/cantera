function setTemperature(d, t)
% SETTEMPERATURE  Set the temperature.
% d = setTemperature(d, t)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :param t:
%     Temperature to be set. Units: K
%

domain_methods(d.dom_id, 61, t);
