function setInitialVolume(r, v0)
% SETINITIALVOLUME  Set the initial reactor volume.
% setInitialVolume(r, v0)
% :param r:
%     Instance of class :mat:func:`Reactor`
% :param v0:
%     Initial volume. Units: m**3
%

reactormethods(4, reactor_hndl(r), v0);
