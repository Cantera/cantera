function mdot = massFlux(d)
% MASSFLUX  Get the mass flux.
% mdot = massFlux(d)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :return:
%     The mass flux in the domain.
%

mdot = domain_methods(d.dom_id, 17);
