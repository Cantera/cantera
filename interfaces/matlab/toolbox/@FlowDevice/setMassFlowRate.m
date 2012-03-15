function setMassFlowRate(f, mdot)
% SETMASSFLOWRATE -
%
if f.type == 1
    k = flowdevicemethods(3, f.index, mdot);
    if k < 0
        error(geterr);
    end
else
    error('Mass flow rate can only be set for mass flow controllers')
end
