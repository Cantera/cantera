function v = volume(r)
% VOLUME  Get the volume of the reactor.
% v = volume(r)
% :param r:
%     Instance of class :mat:func:`Reactor`
% :return:
%     The volume of the reactor contents at the
%     end of the last call to :mat:func:`advance` or :mat:func:`step`.
%     Units: m**3
%

v = reactormethods(24, reactor_hndl(r));
