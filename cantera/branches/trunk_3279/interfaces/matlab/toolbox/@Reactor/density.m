function rho = density(r)
% DENSITY  Get the density of the reactor.
% rho = density(r)
% :param r:
%     Instance of class :mat:func:`Reactor`
% :return:
%     Density of the phase in the input. Units: kg/m**3
%

rho = reactormethods(25, reactor_hndl(r));
