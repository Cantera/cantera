classdef MassFlowController < FlowDevice
    % Create a mass flow controller. ::
    %
    %     >> m = MassFlowController(upstream, downstream, name)
    %
    % Creates an instance of class :mat:class:`FlowDevice` configured to
    % simulate a mass flow controller that maintains a constant mass flow
    % rate independent of upstream or downstream conditions.
    %
    % see also: :mat:class:`FlowDevice`, :mat:class:`Valve`
    %
    % :param upstream:
    %     Upstream :mat:class:`ReactorBase`.
    % :param downstream:
    %     Downstream :mat:class:`ReactorBase`.
    % :param name:
    %     Flow device name (optional; default is ``(none)``).
    % :return:
    %     Instance of class :mat:class:`FlowDevice`.

    methods

        function m = MassFlowController(upstream, downstream, name)
            % Constructor

            arguments
                upstream {mustBeA(upstream, 'ReactorBase')}
                downstream {mustBeA(downstream, 'ReactorBase')}
                name (1,1) string = "(none)"
            end

            m@FlowDevice('MassFlowController', upstream, downstream, name);
        end

    end

end
