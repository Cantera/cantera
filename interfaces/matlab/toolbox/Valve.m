function v = Valve(upstream, downstream)
% VALVE  Create a valve.
% v = Valve(upstream, downstream)
% Create an instance of class :mat:func:`FlowDevice` configured to
% simulate a valve that produces a flow rate proportional to the
% pressure difference between the upstream and downstream reactors.
%
% The mass flow rate [kg/s] is computed from the expression
%
% .. math:: \dot{m} = K(P_{upstream} - P_{downstream})
%
% as long as this produces a positive value.  If this expression is
% negative, zero is returned. Therefore, the :mat:func:`Valve` object
% acts as a check valve - flow is always from the upstream reactor to
% the downstream one.  Note: as currently implemented, the Valve object
% does not model real valve characteristics - in particular, it
% does not model choked flow. The mass flow rate is always assumed
% to be linearly proportional to the pressure difference, no matter how
% large the pressure difference. THIS MAY CHANGE IN A FUTURE
% RELEASE.
%
% see also: :mat:func:`FlowDevice`, :mat:func:`MassFlowController`
%
% :param upstream:
%     Upstream reactor or reservoir
% :param downstream:
%     Downstream Reactor or reservoir
%
% :return:
%     Instance of class :mat:func:`FlowDevice`
%

v = FlowDevice('Valve');
if nargin == 2
    install(v, upstream, downstream)
end
