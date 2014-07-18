function t = temperature(r)
% TEMPERATURE  Get the temperature of the reactor.
% t = temperature(r)
% :param r:
%     Instance of class :mat:func:`Reactor`
% :return:
%     The temperature of the reactor contents at the
%     end of the last call to :mat:func:`advance` or :mat:func:`step`.
%     Units: K
%

t = reactormethods(26, reactor_hndl(r));
