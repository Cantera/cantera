function i = hndl(p)
% HNDL  Get the integer used to access the kernel object.
% i = hndl(p)
% Deprecated in favor of :mat:func:`thermo_hndl`.
%

warning('This function is deprecated in favor of thermo_hndl.m')
i = p.tp_id;
