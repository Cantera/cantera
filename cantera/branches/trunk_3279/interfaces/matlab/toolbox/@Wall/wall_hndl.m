function i = wall_hndl(w)
% WALL_HNDL  Get the integer used to access the kernel object.
% i = wall_hndl(w)
% :param w:
%     Instance of class :mat:func:`Wall`
%     for which the handle is desired.
% :return:
%     Integer used to access the kernel object
%

i = w.index;
