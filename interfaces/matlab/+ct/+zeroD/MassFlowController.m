classdef MassFlowController < ct.zeroD.FlowDevice
    % Create a mass flow controller. ::
    %
    %     >> m = ct.zeroD.MassFlowController(upstream, downstream, name)
    %
    % Creates an instance of class :mat:class:`ct.zeroD.FlowDevice` configured to
    % simulate a mass flow controller that maintains a constant mass flow
    % rate independent of upstream or downstream conditions.
    %
    % see also: :mat:class:`ct.zeroD.FlowDevice`, :mat:class:`ct.zeroD.Valve`
    %
    % :param upstream:
    %     Upstream :mat:class:`ct.zeroD.ReactorBase`.
    % :param downstream:
    %     Downstream :mat:class:`ct.zeroD.ReactorBase`.
    % :param name:
    %     Flow device name (optional; default is ``(none)``).

    methods

        function obj = MassFlowController(upstream, downstream, name)
            arguments
                upstream (1,1) ct.zeroD.ReactorBase
                downstream (1,1) ct.zeroD.ReactorBase
                name (1,1) string = "(none)"
            end

            obj@ct.zeroD.FlowDevice('MassFlowController', upstream, downstream, name);
        end

    end

end
