classdef MassFlowController < FlowDevice
    % Create a mass flow controller. ::
    %
    %     >> m = MassFlowController(upstream, downstream)
    %
    % Creates an instance of class :mat:class:`FlowDevice` configured to
    % simulate a mass flow controller that maintains a constant mass flow
    % rate independent of upstream or downstream conditions. If two reactor
    % objects are supplied as arguments, the controller is installed
    % between the two reactors. Otherwise, the :mat:func:`install` method
    % should be used to install the :mat:class:`MassFlowController` between
    % reactors.
    %
    % see also: :mat:class:`FlowDevice`, :mat:class:`Valve`
    %
    % :param upstream:
    %     Upstream :mat:class:`Reactor` or :mat:class:`Reservoir`
    % :param downstream:
    %     Downstream :mat:class:`Reactor` or :mat:class:`Reservoir`
    % :return:
    %     Instance of class :mat:class:`FlowDevice`
    %

    methods

        function m = MassFlowController(upstream, downstream)
            % Constructor

            m@FlowDevice('MassFlowController');

            if nargin == 2
                m.install(upstream, downstream)
            end

        end

    end

end
