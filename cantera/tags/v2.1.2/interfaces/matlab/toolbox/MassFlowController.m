function m = MassFlowController(upstream, downstream)
%
% MASSFLOWCONTROLLER - Create a mass flow controller connecting two
% reactors / reservoirs.
%
%    m = MassFlowController(upstream, downstream)
%
%    creates an instance of class FlowDevice configured to simulate a
%    mass flow controller that maintains a constant mass flow rate
%    independent of upstream or downstream conditions. If two reactor
%    objects are supplied as arguments, the controller is installed
%    between the two reactors.
%
%    see also: FlowDevice
%
m = FlowDevice(1);
if nargin == 2
    install(m, upstream, downstream)
end
