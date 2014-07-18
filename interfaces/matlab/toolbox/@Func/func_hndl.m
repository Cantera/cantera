function i = func_hndl(f)
% FUNC_HNDL  Get the integer used to access the kernel object.
% i = func_hndl(f)
% :param f:
%     Instance of class :mat:func:`Func`
% :return:
%     The handle of the input function
%

i = f.index;
