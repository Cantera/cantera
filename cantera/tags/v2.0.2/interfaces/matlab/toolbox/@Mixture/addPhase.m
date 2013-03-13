function addPhase(self, phase, moles)
% ADDPHASE - Add a phase to the mixture.
%
%            carbon = importPhase('graphite.cti');
%            addPhase(mix, carbon, 1.0);
%
if ~isa(phase,'ThermoPhase')
    error('phase object of wrong type.');
end
if ~isa(moles,'numeric')
    error('number of moles must be numeric.');
end
if moles < 0.0
    error('negative moles!');
end

iphase = thermo_hndl(phase);
iok = mixturemethods(4, mix_hndl(self), iphase, moles);
if iok < 0
    error('error adding phase');
end
