function v = Valve(upstream, downstream)
    % Create a valve.
    % Creates an instance of class 'FlowDevice' configured to simulate a
    % valve that produces a flow rate proportional to the pressure
    % difference between the upstream or downstream reactors.
    %
    % The mass flow rate [kg/s] is computed from the expression:
    %
    % mdot = K(P_upstream - P_downstream)
    %
    % as long as this produces a positive value. If this expression is
    % negative, zero is returned. Therefore, the 'Valve' object acts as a
    % check valve - flow is always from upstream to downstream.
    %
    % Note: as currently implemented, the valve object does not model real
    % valve characteristics - inparticular, it does not model choked flow.
    % The mass flow rate is always assumed to be linearly proportional to
    % the pressure difference, no matter how large. This MAY change in a
    % future release.
    %
    % :param upstream:
    %    Upstream 'Reactor' or 'Reservoir'.
    % :param downstream:
    %    Downstream 'Reactor' or 'Reservoir.
    % :return:
    %    Instance of class 'FlowDevice'.

    v = FlowDevice('Valve');
    if nargin == 2
        v.install(upstream, downstream)
    end
end
