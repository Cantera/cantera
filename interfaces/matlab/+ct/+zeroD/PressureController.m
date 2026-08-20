classdef PressureController < ct.zeroD.FlowDevice
    % Create a pressure controller. ::
    %
    %     >> p = ct.zeroD.PressureController(upstream, downstream, name)
    %
    % Creates an instance of class :mat:class:`ct.zeroD.FlowDevice` configured to
    % simulate a pressure controller. A :mat:class:`ct.zeroD.PressureController` is
    % designed to be used in conjunction with another primary flow device, typically a
    % :mat:class:`ct.zeroD.MassFlowController`. The primary flow device is installed on
    % the inlet of the reactor, and the corresponding
    % :mat:class:`ct.zeroD.PressureController` is installed on the outlet of the
    % reactor. The mass flow rate through the :mat:class:`ct.zeroD.PressureController`
    % is equal to the mass flow rate of the primary device, plus a correction that
    % depends on the pressure difference:
    %
    % .. math:: \dot{m} = \dot{m}_{primary} + K_v(P_{upstream} - P_{downstream})
    %
    % as long as this produces a positive value. If this expression is negative, zero
    % is returned. The primary device is set using
    % :mat:meth:`ct.zeroD.FlowDevice.setPrimary` and must be set before the mass flow
    % rate can be evaluated.
    %
    % see also: :mat:class:`ct.zeroD.FlowDevice`,
    % :mat:class:`ct.zeroD.MassFlowController`, :mat:class:`ct.zeroD.Valve`
    %
    % :param upstream:
    %     Upstream :mat:class:`ct.zeroD.ReactorBase`.
    % :param downstream:
    %     Downstream :mat:class:`ct.zeroD.ReactorBase`.
    % :param name:
    %     Flow device name (optional; default is ``(none)``).

    properties (SetAccess = public)

        % Proportionality constant :math:`K_v` in kg/s/Pa between the pressure drop
        % and the mass flow rate correction added to the primary device's mass flow
        % rate.
        pressureCoeff

    end

    methods

        function obj = PressureController(upstream, downstream, name)
            arguments
                upstream (1,1) ct.zeroD.ReactorBase
                downstream (1,1) ct.zeroD.ReactorBase
                name (1,1) string = "(none)"
            end

            obj@ct.zeroD.FlowDevice('PressureController', upstream, downstream, name);
        end

        %% PressureController Get Methods

        function k = get.pressureCoeff(obj)
            k = obj.deviceCoefficient;
        end

        %% PressureController Set Methods

        function set.pressureCoeff(obj, k)
            obj.deviceCoefficient = k;
        end

    end

end
