function n = hndl(a)
% HNDL  Get the integer used to access the kernel object.
% n = hndl(a)
% Integer used to access kernel objects. Deprecated in favor of
% :mat:func:`trans_hndl`.
%

warning('This function is deprecated in favor of trans_hndl')
n = a.id;
