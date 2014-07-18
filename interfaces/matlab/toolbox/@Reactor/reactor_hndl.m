function i = reactor_hndl(r)
% REACTOR_HNDL  Get the integer used to access the kernel object.
% i = reactor_hndl(r)
% :param r:
%     Instance of class :mat:func:`Reactor`
%     for which the handle is desired.
% :return:
%     Integer used to access the kernel object
%

i = r.index;
