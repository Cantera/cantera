function setMassFlowRate(f, mdot)
% SETMASSFLOWRATE  Set the mass flow rate to a constant value.
% setMassFlowRate(f, mdot)
%
% See also: :mat:func:`MassFlowController`
%
% :param f:
%     Instance of class :mat:func:`MassFlowController`
% :param mdot:
%     Mass flow rate
%
if strcmp(f.type, 'MassFlowController')
    k = flowdevicemethods(10, f.index, mdot);
    if k < 0
        error(geterr);
    end
else
    error('Mass flow rate can only be set for mass flow controllers')
end
