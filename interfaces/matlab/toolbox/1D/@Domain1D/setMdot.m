function setMdot(d, mdot)
% SETMDOT  Set the mass flow rate.
% d = setMdot(d, mdot)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :param mdot:
%     Mass flow rate
%

domain_methods(d.dom_id, 60, mdot);
