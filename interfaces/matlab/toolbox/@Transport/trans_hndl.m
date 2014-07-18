function n = trans_hndl(a)
% TRANS_HNDL  Get the integer used to access the kernel object.
% n = trans_hndl(a)
% :param a:
%     Instance of class :mat:func:`Transport` (or another
%     object derived from Transport)
%     for which the handle is desired.
% :return:
%     Integer used to access the kernel object
%

n = a.id;
