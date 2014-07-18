function addPhase(self, phase, moles)
% ADDPHASE - Add a phase to the mixture.
%
%            carbon = importPhase('graphite.cti');
%            addPhase(mix, carbon, 1.0);
%
if ~isa(phase,'ThermoPhase')
    error('Phase object of wrong type.');
end
if ~isa(moles,'numeric')
    error('Number of moles must be numeric.');
end
if moles < 0.0
    error('Negative moles!');
end

iphase = thermo_hndl(phase);
iok = mixturemethods(4, mix_hndl(self), iphase, moles);
if iok < 0
    error('Error adding phase');
end
