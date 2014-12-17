function p = pressure(r)
% PRESSURE  Get the pressure of the reactor.
% p = pressure(r)
% :param r:
%     Instance of class :mat:func:`Reactor`
% :return:
%     The pressure of the reactor contents at the
%     end of the last call to :mat:func:`advance` or :mat:func:`step`.
%     Units: Pa
%

p = reactormethods(29, reactor_hndl(r));
