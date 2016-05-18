function install(s, reactor)
% INSTALL  Install a ReactorSurface in a Reactor.
% install(s, reactor)
% :param s:
%     Instance of class :mat:func:`ReactorSurface`
% :param reactor:
%     Instance of class :mat:func:`Reactor`
%

s.reactor = reactor;
reactorsurfacemethods(4, reactorsurface_hndl(s), reactor_hndl(reactor));
