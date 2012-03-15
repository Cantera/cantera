function m = Valve(upstream, downstream)
%
% VALVE - Create a valve connecting two reactors / reservoirs.
%
%    m = Valve(upstream, downstream)
%
%    creates an instance of class FlowDevice configured to simulate a
%    valve that produces a flow rate proportional to the pressure
%    difference between the uupstream and downstream reactors. If two reactor
%    objects are supplied as arguments, the valve is installed
%    between the two reactors.
%
%    The mass flow rate [kg/s] is computed from the expression
%
%        mdot = K ( P_upstream - P_downstream )
%
%    as long as this produces a positive value.  If this expression is
%    negative, zero is returned. Therefore, the Valve object acts as a
%    check valve - flow is always from the upstream reactor to the
%    downstream one.  Note: as currently implemented, the Valve object
%    does not model real valve characteristics - in particular, it
%    does not model choked flow. The mass flow rate is always assumed
%    to be linearly proportional to the mass flow rate, no matter how
%    large the pressure difference. THIS MAY CHANGE IN A FUTURE
%    RELEASE.
%
%    see also: FlowDevice, MassFlowController
%
m = FlowDevice(3);
if nargin == 2
    install(m, upstream, downstream)
end
