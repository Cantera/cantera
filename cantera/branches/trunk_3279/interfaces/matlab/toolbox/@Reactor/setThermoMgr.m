function setThermoMgr(r, t)
% SETTHERMOMGR  Set the thermodynamics manager.
% setThermoMgr(r, t)
% This method is used internally during Reactor initialization, but
% is usually not called by users.
%
% :param r:
%     Instance of class :mat:func:`Reactor`
% :param t:
%     Instance of class :mat:func:`ThermoPhase`, or another object
%     containing an instance of that class.
%

if ~isa(t,'ThermoPhase')
    error('wrong object type');
end

reactormethods(6, reactor_hndl(r), thermo_hndl(t));
