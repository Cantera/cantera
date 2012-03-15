function setThermoMgr(r, t)
% SETTHERMOMGR - set the thermo manager
%
if ~isa(t,'ThermoPhase')
    error('wrong object type');
end

reactormethods(6, reactor_hndl(r), thermo_hndl(t));
