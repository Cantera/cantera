function i = thermo_hndl(tp)
% THERMO_HNDL  Get the integer used to access the kernel object.
% i = thermo_hndl(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object deriving from ThermoPhase)
%     for which the handle is desired.
% :return:
%     Integer used to access the kernel object.
%

i = tp.tp_id;
