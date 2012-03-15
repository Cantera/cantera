function mdot = massFlowRate(f, time)
% MASSFLOWRATE - mass flow rate in kg/s
%
mdot = flowdevicemethods(21, f.index, time);
