function addPhase(self, phase, moles)
% ADDPHASE  Add a phase to a mixture.
% addPhase(self, phase, moles)
% :param self:
%     Instance of class :mat:func:`Mixture` to which phases should be
%     added
% :param phase:
%     Instance of class :mat:func:`ThermoPhase` which should be added
% :param moles:
%     Number of moles of the ``phase`` to be added to this mixture.
%     Units: kmol
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
