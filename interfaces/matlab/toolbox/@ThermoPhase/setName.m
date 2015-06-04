function setName(tp, name)
% SETNAME  Set the name of the phase.
% setName(tp, name)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param name:
%     String, name of the phase
%

if isa(name,'char')
    phase_set(thermo_hndl(tp), 32, name);
else
    error('name must be a string.');
end
