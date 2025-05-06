classdef MassFlowController < FlowDevice
    % Create a mass flow controller. ::
    %
    %     >> m = MassFlowController(upstream, downstream, name)
    %
    % Creates an instance of class :mat:class:`FlowDevice` configured to
    % simulate a mass flow controller that maintains a constant mass flow
    % rate independent of upstream or downstream conditions. If two reactor
    % objects are supplied as arguments, the controller is installed
    % between the two reactors.
    %
    % see also: :mat:class:`FlowDevice`, :mat:class:`Valve`
    %
    % :param upstream:
    %     Upstream :mat:class:`Reactor` or :mat:class:`Reservoir`.
    % :param downstream:
    %     Downstream :mat:class:`Reactor` or :mat:class:`Reservoir`.
    % :param name:
    %     Flow device name (optional; default is ``(none)``).
    % :return:
    %     Instance of class :mat:class:`FlowDevice`.

    methods

        function m = MassFlowController(upstream, downstream, name)
            % Constructor

            if nargin < 3
                name = '(none)';
            end

            m@FlowDevice('MassFlowController', upstream, downstream, name);
        end

    end

end
