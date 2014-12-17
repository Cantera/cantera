function i = hndl(k)
% HNDL  Get the integer used to access kernel objects.
% i = hndl(k)
% Integer used to access kernel objects. Deprecated in favor of
% :mat:func:`kinetics_hndl`.
%

warning('This function is deprecated in favor of kinetics_hndl.m')
i = k.id;
