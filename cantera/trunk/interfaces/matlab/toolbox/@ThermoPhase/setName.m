function setName(tp, name)
% SETNAME - Set the name of the phase.
%

if isa(name,'char')
    phase_set(thermo_hndl(tp), 32, name);
else
    error('name must be a string.');
end
