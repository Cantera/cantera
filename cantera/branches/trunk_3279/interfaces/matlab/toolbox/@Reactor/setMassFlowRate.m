function setMassFlowRate(r, mdot)
% SETMASSFLOWRATE  Set the mass flow rate.
% setMassFlowRate(r, mdot)
% :param r:
%     Instance of class :mat:func:`Reactor`
% :param mdot:
%     Mass flow rate. Units: kg/s
%

reactormethods(10, reactor_hndl(r), mdot);
