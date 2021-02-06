function mdot = massFlowRate(f, time)
% MASSFLOWRATE  Get the mass flow rate at a given time.
% mdot = massFlowRate(f, time)
% :param f:
%     Instance of class :mat:func:`MassFlowController`
% :param time:
%     Time at which the mass flow rate is desired
% :return:
%     The mass flow rate through the :mat:func:`FlowDevice` at the given time
%

if nargin == 1
    mdot = flowdevicemethods(21, f.index);
else
    warning(['"time" argument to massFlowRate is deprecated and will be' ...
            ' removed after Cantera 2.5.'])
    mdot = flowdevicemethods(21, f.index, time);
end
