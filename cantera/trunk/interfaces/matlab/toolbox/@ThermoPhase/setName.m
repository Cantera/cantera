function setName(self, name)
% SETNAME - Set the name of the phase.
%
if isa(name,'char')
    phase_set(thermo_hndl(self), 32, name);
else
    error('name must be a string.');
end
