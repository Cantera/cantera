function nm = name(tp)
% NAME  Get the name of the phase.
% nm = name(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :return:
%     String name of the input phase
%

nm = phase_get(thermo_hndl(tp), 42);
