classdef MassFlowController < ct.FlowDevice
    % Create a mass flow controller. ::
    %
    %     >> m = ct.MassFlowController(upstream, downstream, name)
    %
    % Creates an instance of class :mat:class:`ct.FlowDevice` configured to
    % simulate a mass flow controller that maintains a constant mass flow
    % rate independent of upstream or downstream conditions.
    %
    % see also: :mat:class:`ct.FlowDevice`, :mat:class:`ct.Valve`
    %
    % :param upstream:
    %     Upstream :mat:class:`ct.ReactorBase`.
    % :param downstream:
    %     Downstream :mat:class:`ct.ReactorBase`.
    % :param name:
    %     Flow device name (optional; default is ``(none)``).

    methods

        function m = MassFlowController(upstream, downstream, name)
            arguments
                upstream (1,1) ct.ReactorBase
                downstream (1,1) ct.ReactorBase
                name (1,1) string = "(none)"
            end

            m@ct.FlowDevice('MassFlowController', upstream, downstream, name);
        end

    end

end
