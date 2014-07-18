function m = mass(r)
% MASS  Get the mass of the reactor.
% m = mass(r)
% :param r:
%     Instance of class :mat:func:`Reactor`
% :return:
%     The mass of the reactor contents at the
%     end of the last call to :mat:func:`advance` or :mat:func:`step`.
%     The mass is retrieved from the solution vector. Units: kg
%

m = reactormethods(23, reactor_hndl(r));
