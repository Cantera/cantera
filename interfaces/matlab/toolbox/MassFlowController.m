function m = MassFlowController(upstream, downstream)
% MASSFLOWCONTROLLER  Create a mass flow controller.
% m = MassFlowController(upstream, downstream)
% Creates an instance of class :mat:func:`FlowDevice` configured to
% simulate a mass flow controller that maintains a constant mass flow
% rate independent of upstream or downstream conditions. If two reactor
% objects are supplied as arguments, the controller is installed
% between the two reactors. Otherwise, the :mat:func:`install` method
% should be used to install the :mat:func:`MassFlowController` between
% reactors.
%
% see also: :mat:func:`FlowDevice`, :mat:func:`Valve`
%
% :param upstream:
%     Upstream :mat:func:`Reactor` or :mat:func:`Reservoir`
% :param downstream:
%     Downstream :mat:func:`Reactor` or :mat:func:`Reservoir`
% :return:
%     Instance of class :mat:func:`FlowDevice`
%

m = FlowDevice('MassFlowController');
if nargin == 2
    install(m, upstream, downstream)
end
