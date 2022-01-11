function m = MassFlowController(upstream, downstream)
    % Create a mass flow controller.
    % Creates an instance of class 'FlowDevice' configured to simulate a
    % mass flow controller that maintains a constant mass flow rate
    % independent of upstream or downstream conditions. If two reactor
    % objects are supplied as arguments, the controller is installed
    % between the two reactors. Otherwise, the 'install' method should be
    % used to install the mass flow controller between reactors.
    %
    % :param upstream:
    %    Upstream 'Reactor' or 'Reservoir'.
    % :param downstream:
    %    Downstream 'Reactor' or 'Reservoir.
    % :return:
    %    Instance of class 'FlowDevice'.

    m = FlowDevice('MassFlowController');
    if nargin == 2
        m.install(upstream, downstream)
    end
end
