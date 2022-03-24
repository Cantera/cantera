function mdot = massFlowRate(f)
% MASSFLOWRATE  Get the mass flow rate.
% mdot = massFlowRate(f)
% :param f:
%     Instance of class :mat:func:`MassFlowController`
% :return:
%     The mass flow rate through the :mat:func:`FlowDevice` at the current time
%

mdot = flowdevicemethods(21, f.index);
