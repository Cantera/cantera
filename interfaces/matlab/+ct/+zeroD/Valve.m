classdef Valve < ct.zeroD.FlowDevice
    % Create a valve. ::
    %
    %     >> v = ct.zeroD.Valve(upstream, downstream, name)
    %
    % Create an instance of class :mat:class:`ct.zeroD.FlowDevice` configured to
    % simulate a valve that produces a flow rate proportional to the
    % pressure difference between the upstream and downstream reactors.
    %
    % The mass flow rate [kg/s] is computed from the expression
    %
    % .. math:: \dot{m} = K(P_{upstream} - P_{downstream})
    %
    % as long as this produces a positive value.  If this expression is
    % negative, zero is returned. Therefore, the :mat:class:`ct.zeroD.Valve` object
    % acts as a check valve - flow is always from the upstream reactor to
    % the downstream one.  Note: as currently implemented, the Valve object
    % does not model real valve characteristics - in particular, it
    % does not model choked flow. The mass flow rate is always assumed
    % to be linearly proportional to the pressure difference, no matter how
    % large the pressure difference. THIS MAY CHANGE IN A FUTURE RELEASE.
    %
    % see also: :mat:class:`ct.zeroD.FlowDevice`, :mat:class:`ct.zeroD.MassFlowController`
    %
    % :param upstream:
    %     Upstream :mat:class:`ct.zeroD.ReactorBase`.
    % :param downstream:
    %     Downstream :mat:class:`ct.zeroD.ReactorBase`.
    % :param name:
    %     Flow device name (optional; default is ``(none)``).

    methods

        function obj = Valve(upstream, downstream, name)
            arguments
                upstream (1,1) ct.zeroD.ReactorBase
                downstream (1,1) ct.zeroD.ReactorBase
                name (1,1) string = "(none)"
            end

            obj@ct.zeroD.FlowDevice('Valve', upstream, downstream, name);
        end

    end

end
