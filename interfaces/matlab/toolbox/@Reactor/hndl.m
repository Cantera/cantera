function i = hndl(r)
% HNDL  Get the integer used to access the kernel object.
% i = hndl(r)
% Integer used to access kernel objects. Deprecated in favor of
% :mat:func:`reactor_hndl`.
%

warning('This function is deprecated in favor of reactor_hndl.m')
i = r.index;
